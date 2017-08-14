#pragma once

#include <iostream>
#include <vector>
#include <assert.h>
#include <utility>
#include <chrono>
#include <ctime>
#include <random>
#include <vector>
#include "helper_functions.h"

//! Evolutionary Algorithms
namespace ea
{
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
	//! Random number generator
	std::mt19937 generator(std::random_device{}());

	//! Template Class for Solvers
	template<template<typename> class S, typename T, typename F, typename C> class Solver;

	//! Base Class for Evolutionary Algorithms
	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	class Solver_base
	{
	public:
		//! Displays the results
		void display_results();
		//! Solve wrapper function for Solvers, used for benchmarks
		std::vector<T> solver_bench()
		{
			std::cout << static_cast<Derived*>(this)->type << " used as solver." << "\n";
			if (solver_struct.tol > std::abs(f(min_cost)))
			{
				timer = 0;
				display_results();
				return min_cost;
			}
			else
			{
				//! Time the computation
				const std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
				static_cast<Derived*>(this)->run_algo();
				const std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
				const std::chrono::duration<double> elapsed_seconds = end - start;
				timer = elapsed_seconds.count();
				display_results();
				//! Return minimum cost individual
				return min_cost;
			}
		}
	protected:
		//! Constructor
		Solver_base(const S<T>& i_solver_struct, const F& i_f, const C& i_c)
			: solver_struct{ i_solver_struct },
			individuals{ init_individuals() },
			min_cost{ individuals[0] }, iter{ 0 }, f{ i_f }, c{ i_c }, solved_flag{ false }, timer{ 0 }

		{
            distribution = std::uniform_real_distribution<T>(0.0,1.0);
			find_min_cost();
		}
		//! Internal copy of the structure used for parameters of the algorithm
		const S<T> solver_struct;
		/*! /brief The fitness and constraints functions are copied
		so that even if they are not available in the current scope, the solver will still execute properly.
		At the same time, std::function could be have used, thus eliminating the need for template parameters F and C.
		However, that comes at a runtime cost, since calls to the functions would be virtual and there is a possibility
		that allocation could happen on the heap.*/
		//! Population
		std::vector<std::vector<T>> individuals;
		//! Best solution / lowest fitness
		std::vector<T> min_cost;
		//! Last iteration to solution
		size_t iter;
		//! Copy of the fitness function passed as a lambda
		F f;
		//! Copy of the constraints function passed as a lambda
		C c;
		//! Uniform real distribution
		std::uniform_real_distribution<T> distribution;
		//! A flag which determines if the solver has already solved the problem
		bool solved_flag;
		//! Timer
		T timer;
		//! Returns a randomised individual using the initial decision variables and standard deviation
		std::vector<T> randomise_individual();
		//! Initialises the population by randomising aroung the decision variables using the given standard deviation
		std::vector<std::vector<T>> init_individuals();
		//! Find the minimum cost individual of the fitness function for the population
		void find_min_cost();
	};

	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	std::vector<T> Solver_base<Derived, S, T, F, C>::randomise_individual()
	{
		std::vector<T> individual = solver_struct.decision_variables;
		T epsilon = 0;
		for (size_t j = 0; j < solver_struct.ndv; ++j)
		{
			std::normal_distribution<T> ndistribution(0, solver_struct.stdev[j]);
			epsilon = ndistribution(generator);
			individual[j] = individual[j] + epsilon;
		}
		return individual;
	}

	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	std::vector<std::vector<T>> Solver_base<Derived, S, T, F, C>::init_individuals()
	{
		generator.discard(700000);
		std::vector<std::vector<T>> individuals(solver_struct.npop, std::vector<T>(solver_struct.ndv));
		for (auto& p : individuals)
		{
			p = randomise_individual();
			//! Check population constraints
			while (!c(p))
			{
				p = randomise_individual();
			}
		}
		return individuals;
	}

	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	void Solver_base<Derived, S, T, F, C>::find_min_cost()
	{
		for (const auto& p : individuals)
		{
			if (f(min_cost) > f(p))
			{
                min_cost = p;
			}
		}
	}

	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	void Solver_base<Derived, S, T, F, C>::display_results()
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

	/*!
	/brief Solver wrapper function, interface to solvers : free function used for benchmarks
	*/
	template<typename F, typename C, template<typename> class S, typename T>
	std::vector<T> solve(const F& f, const C& c, const S<T>& solver_struct)
	{
		Solver<S, T, F, C> solver{ solver_struct, f, c };
		return solver.solver_bench();
	}
}