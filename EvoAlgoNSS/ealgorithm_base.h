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
	//! The floating point number type used for type deduction
	using fp_type = T;
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
protected:
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
};

//! Template Class for Solvers
template<typename S, typename F, typename C> class Solver;

//! Base Class for Evolutionary Algorithms
template<typename T, typename F, typename C>
class Solver_base
{
public:
	//! Displays the results
	void display_results();
	//! solve wrapper function for Solvers, used for benchmarks
	template<typename S> std::vector<T> solver_bench(Solver<S, F, C>& solver, const S& solver_struct)
	{
		std::cout << solver.type << " used as solver." << "\n";
		if (solver_struct.tol > std::abs(f(solver.min_cost)))
		{
			timer = 0;
			solver.display_results();
			return solver.min_cost;
		}
		else
		{
			//! Time the computation
			std::chrono::time_point<std::chrono::system_clock> start, end;
			start = std::chrono::system_clock::now();
			solver.run_algo();
			end = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = end - start;
			std::time_t end_time = std::chrono::system_clock::to_time_t(end);
			timer = elapsed_seconds.count();
			solver.display_results();
			//! Return minimum cost individual
			return solver.min_cost;
		}
	}
protected:
	//! Constructor
	Solver_base(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev, const T& tol, const F& i_f, const C& i_c)
		: individuals{ init_individuals(decision_variables, npop, stdev) }, min_cost{ individuals[0] }, iter{ 0 }, f{ i_f }, c{ i_c }, solved_flag{ false }, timer{ 0 }
	{
		distribution = std::uniform_real_distribution<T>::uniform_real_distribution(0.0, 1.0);
		population_constraints_checker(decision_variables, stdev);
		find_min_cost();
	}
	/*! The fitness and constraints functions are copied 
	so that even if they are not available in the current scope, the solver will still execute properly.
	At the same time, std::function could be have used, thus eliminating the need for template parameters F and C.
	However, that comes at a runtime cost, since calls to the functions would be virtual and there is a possibility
	that allocation could happen on the heap.!*/
	//! Copy of the fitness function passed as a lambda
	F f;
	//! Copy of the constraints function passed as a lambda
	C c;
	//! Population
	std::vector<std::vector<T>> individuals;
	//! Last iteration to solution
	size_t iter;
	//! Best solution / lowest fitness
	std::vector<T> min_cost;
	//! Random number generator
	std::random_device generator;
	//! Uniform real distribution
	std::uniform_real_distribution<T> distribution;
	//! A flag which determines if the solver has already solved the problem
	bool solved_flag;
	//! Timer
	T timer;
	//! Returns a randomised individual using the initial decision variables and standard deviation
	std::vector<T> randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev);
	//! Initialises the population by randomising aroung the decision variables using the given standard deviation
	std::vector<std::vector<T>> init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev);
	//! Check population constraints
	void population_constraints_checker(const std::vector<T>& decision_variables, const std::vector<T>& stdev);
	//! Find the minimum cost individual of the fitness function for the population
	void find_min_cost();
};

template<typename T, typename F, typename C>
std::vector<T> Solver_base<T, F, C>::randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev)
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

template<typename T, typename F, typename C>
std::vector<std::vector<T>> Solver_base<T, F, C>::init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev)
{
	const auto& ndv = decision_variables.size();
	std::vector<std::vector<T>> individuals(npop, std::vector<T>(ndv));
	for (auto& p : individuals)
	{
		p = randomise_individual(decision_variables, stdev);
	}
	return individuals;
}

template<typename T, typename F, typename C>
void Solver_base<T, F, C>::population_constraints_checker(const std::vector<T>& decision_variables, const std::vector<T>& stdev)
{
	for (auto& p : individuals)
	{
		while (!c(p))
		{
			p = randomise_individual(decision_variables, stdev);
		}
	}
};

template<typename T, typename F, typename C>
void Solver_base<T, F, C>::find_min_cost()
{
	for (const auto& p : individuals)
	{
		if (f(min_cost) > f(p))
		{
			min_cost = p;
		}
	}
}

template<typename T, typename F, typename C>
void Solver_base<T, F, C>::display_results()
{
	if (!solved_flag)
	{
		std::cout << "!!!The optimisation problem was not solved.!!!" << "\n";
	}
	else
	{
		std::cout << "!!!Problem was successfully solved.!!!" << "\n";
	}
	std::cout << "Optimum solution: " << min_cost << " Fitness Value: " << f(min_cost) << "\n";
	std::cout << "Population: " << individuals.size() << " Solved at iteration: " << iter << "\n";
	std::cout << "Elapsed time in seconds: " << timer << "\n";
}

//! Solver wrapper function, interface to solvers : free function used for benchmarks
template<typename F, typename C, typename S, typename T = S::fp_type>
std::vector<T> solve(const F& f, const C& c, const S& solver_struct)
{
	Solver<S, F, C> solver(solver_struct, f, c);
	return solver.solver_bench(solver, solver_struct);
}