#pragma once

#include <vector>
#include <assert.h>
#include <utility>
#include <chrono>
#include <ctime>
#include <random>
#include <vector>
#include "helper_functions.h"

//! Evolutionary algorithm stucture base
template<typename T>
struct EA_base
{
public:
	//! Constructor
	EA_base(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
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
	//! The floating point number type used for type deduction
	typedef T fp_type;
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
};

//! Template Class for Solvers
template<typename T, typename S> class Solver;

//! Base Class for Evolutionary Algorithms
template<typename T>
class Solver_base
{
public:
	//! Constructor
	template<typename F, typename C> Solver_base(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev, F f, C c)
		: individuals{ init_individuals(decision_variables, npop, stdev) }, min_cost{ individuals[0] }, iter{ 0 }, solved_flag { false }
	{
		population_constraints_checker(decision_variables, stdev, c);
		find_min_cost(f);
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	//! Returns the minimum cost solution
	std::vector<T> ret_min_cost() { return min_cost; };
	//! Returns the population size
	size_t ret_npop() const { return individuals.size(); };
	//! Returns the fitness cost of the minimum cost solution
	T ret_fitness_cost() const { return fitness_cost; };
	//! Returns the number of iterations that were executed
	size_t ret_iter() const { return iter; };
	//! Displays the results
	void display_results();
protected:
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
	//! A flag which determines if the solver has already solved the problem
	bool solved_flag;
	//! Returns a randomised individual using the initial decision variables and standard deviation
	std::vector<T> randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev);
	//! Initialises the population by randomising aroung the decision variables using the given standard deviation
	std::vector<std::vector<T>> init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev);
	//! Check population constraints
	template<typename C> void population_constraints_checker(const std::vector<T>& decision_variables, const std::vector<T>& stdev, C c);
	//! Find the minimum cost individual of the fitness function for the population
	template<typename F> void find_min_cost(F f);
};

template<typename T>
std::vector<T> Solver_base<T>::randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev)
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

template<typename T>
std::vector<std::vector<T>> Solver_base<T>::init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev)
{
	const auto& ndv = decision_variables.size();
	std::vector<std::vector<T>> individuals(npop, std::vector<T>(ndv));
	for (auto& p : individuals)
	{
		p = randomise_individual(decision_variables, stdev);
	}
	return individuals;
}

template<typename T>
template<typename C>
void Solver_base<T>::population_constraints_checker(const std::vector<T>& decision_variables, const std::vector<T>& stdev, C c)
{
	for (auto& p : individuals)
	{
		while (!c(p))
		{
			p = randomise_individual(decision_variables, stdev);
		}
	}
};

template<typename T>
template<typename F>
void Solver_base<T>::find_min_cost(F f)
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

template<typename T>
void Solver_base<T>::display_results()
{
	std::cout << "Optimum solution: " << min_cost << " Fitness Value: " << fitness_cost << "\n";
	std::cout << "Population: " << individuals.size() << " Solved at iteration: " << iter << "\n";
}

//! solve wrapper function for Solvers, used for benchmarks
template<typename F, typename C, typename S, typename T = S::fp_type>
std::vector<T> solve(F f, C c, const S& solver_struct)
{
	Solver<T, S> solver(solver_struct, f, c);
	std::cout << solver.type << " used as solver" << "\n";
	if (solver_struct.tol > std::abs(solver.ret_fitness_cost()))
	{
		T timer = 0;
		solver.display_results();
		std::cout << "Elapsed time in seconds: " << timer << "\n";
		return solver.ret_min_cost();
	}
	else
	{
		//! Time the computation
		std::chrono::time_point<std::chrono::system_clock> start, end;
		start = std::chrono::system_clock::now();
		solver.run_algo(f, c);
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);
		T timer = elapsed_seconds.count();
		solver.display_results();
		std::cout << "Elapsed time in seconds: " << timer << "\n";
		//! Return minimum cost individual
		return solver.ret_min_cost();
	}
}