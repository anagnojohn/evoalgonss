#pragma once

#include <vector>
#include <assert.h>
#include <utility>
#include <chrono>
#include <ctime>
#include <random>
#include <vector>
#include <tuple>
#include "helper_functions.h"

//! Base Class for Evolutionary Algorithms
template<typename T, typename S>
class Solver_base
{
public:
	//! Type of floating number, used for type deduction
	typedef T value_type;
	//! Displays the results
	void display_results();
	//! Return ndv used for checks of the number of variables used outside of the solver
	size_t ret_ndv() const { return ndv; };
	//! solve wrapper function for Solvers, used for benchmarks
	template<typename F, typename C> std::vector<T> solve(F f, C c);
protected:
	//! Constructor
	Solver_base(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: decision_variables{ i_decision_variables }, stdev{ i_stdev }, npop{ i_npop }, tol{ i_tol }, iter_max{ i_iter_max }, ndv{ i_decision_variables.size() }
	{
		assert(decision_variables.size() > 0);
		assert(decision_variables.size() == stdev.size());
		for (const auto& p : stdev)
		{
			assert(p > 0);
		}
		assert(npop > 0);
		assert(tol > 0);
		assert(iter_max > 0);
		
	}
	//! Decision Variables
	const std::vector<T> decision_variables;
	//! Standard deviation of the decision variables
	const std::vector<T> stdev;
	//! Size of the population
	const size_t npop;
	//! Tolerance
	const T tol;
	//! Number of maximum iterations
	const size_t iter_max;
	//! Number of decision variables
	const size_t ndv;
	//! Population
	std::vector<std::vector<T>> individuals;
	//! Best fitness
	T fitness_cost;
	//! Last iteration to solution
	size_t iter;
	//! Best solution
	std::vector<T> min_cost;
	//! Random number generator
	std::random_device generator;
	//! Uniform real distribution
	std::uniform_real_distribution<T> distribution;
	//! name of the algorithm
	const std::string type = "";
	//! Start time of algorithm
	std::chrono::time_point<std::chrono::system_clock> start;
	//! End time algorithm
	std::chrono::time_point<std::chrono::system_clock> end;
	//! Returns a randomised individual using the initial decision variables and standard deviation
	std::vector<T> randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev);
	//! Initialises the population by randomising aroung the decision variables using the given standard deviation
	std::vector<std::vector<T>> init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev);
	//! Check population constraints
	template<typename C> void population_constraints_checker(const std::vector<T>& decision_variables, const std::vector<T>& stdev, C c);
	//! Find the minimum cost individual of the fitness function for the population
	template<typename F> void find_min_cost(F f);
	template<typename F, typename C> bool initial_solution_checker(F f, C c);
	//! run_algo is the function that is run by each of the algorithms (overloaded in the algorithm classes)
	template<typename F, typename C> void run_algo(F f, C c) {};
};

template<typename T, typename S>
std::vector<T> Solver_base<T, S>::randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev)
{
	std::random_device generator;
	std::vector<T> individual = decision_variables;
	const auto& ndv = decision_variables.size();
	T epsilon;
	for (auto j = 0; j < ndv; ++j)
	{
		std::normal_distribution<T> ndistribution(0, stdev[j]);
		epsilon = std::abs(ndistribution(generator));
		individual[j] = individual[j] + epsilon;
	}
	return individual;
}

template<typename T, typename S>
std::vector<std::vector<T>> Solver_base<T, S>::init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev)
{
	const auto& ndv = decision_variables.size();
	std::vector<std::vector<T>> individuals(npop, std::vector<T>(ndv));
	for (auto& p : individuals)
	{
		p = randomise_individual(decision_variables, stdev);
	}
	return individuals;
}

template<typename T, typename S>
template<typename C>
void Solver_base<T, S>::population_constraints_checker(const std::vector<T>& decision_variables, const std::vector<T>& stdev, C c)
{
	for (auto& p : individuals)
	{
		while (!c(p))
		{
			p = randomise_individual(decision_variables, stdev);
		}
	}
};

template<typename T, typename S>
template<typename F>
void Solver_base<T, S>::find_min_cost(F f)
{
	for (const auto& p : individuals)
	{
		if (f(min_cost) > f(p))
		{
			min_cost = p;
		}
	}
	fitness_cost = f(min_cost);
}

template<typename T, typename S>
void Solver_base<T, S>::display_results()
{
	std::cout << "Optimum solution: " << min_cost << " Fitness Value: " << fitness_cost << "\n";
	std::cout << "Population: " << individuals.size() << " Solved at iteration: " << iter << "\n";
	
}

template<typename T, typename S>
template<typename F, typename C>
bool Solver_base<T, S>::initial_solution_checker(F f, C c)
{
	std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
	distribution = i_distribution;
	individuals = init_individuals(decision_variables, npop, stdev);
	min_cost = decision_variables;
	iter = 0;
	population_constraints_checker(decision_variables, stdev, c);
	find_min_cost(f);
	if (tol > std::abs(fitness_cost))
	{
		T timer = 0;
		display_results();
		std::cout << "Elapsed time in seconds: " << timer << "\n";
		return true;
	}
	else
	{
		return false;
	}
}

template<typename T, typename S>
template<typename F, typename C>
std::vector<T> Solver_base<T, S>::solve(F f, C c)
{
	std::cout << static_cast<S*>(this)->type << " used as solver" << "\n";
	if (initial_solution_checker(f, c))
	{
		return min_cost;
	}
	else
	{
		//! Time the computation
		static_cast<S*>(this)->run_algo(f, c);
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);
		T timer = elapsed_seconds.count();
		display_results();
		std::cout << "Elapsed time in seconds: " << timer << "\n";
		//! Return minimum cost individual
		return min_cost;
	}
}

//! Solver non-member function, creates a Solver and executes it
template<typename F, typename C, typename S, typename T = typename S::value_type>
std::vector<T> solve(F f, C c, S& solver)
{
	auto result = solver.solve(f, c);
	return result;
}
