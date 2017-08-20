#pragma once

#include <iostream>
#include <vector>
#include <assert.h>
#include <utility>
#include <chrono>
#include <ctime>
#include <random>
#include <vector>
#include <fstream>
#include <sstream>
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
		//! Initial Decision Variables
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
		//! Use penalty function or not
		const bool use_penalty_method;
		//! Constraints type
		const Constraints_type constraints_type;
		//! Print to output or not
		const bool print_to_output;
		//! Print to file
		const bool print_to_file;
	protected:
		//! Constructor
		EA_base(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max,
			const bool& i_use_penalty_method, const Constraints_type& i_constraints_type, const bool& i_print_to_output, const bool& i_print_to_file)
			: decision_variables{ i_decision_variables }, stdev{ i_stdev }, npop{ i_npop }, tol{ i_tol }, iter_max{ i_iter_max }, ndv{ i_decision_variables.size() },
			use_penalty_method{ i_use_penalty_method }, constraints_type{ i_constraints_type }, print_to_output{ i_print_to_output }, print_to_file{ i_print_to_file }
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
	
	std::random_device rd;
	//! Random number generator
	std::mt19937_64 generator(rd());

	//! Template Class for Solvers
	template<template<typename> class S, typename T, typename F, typename C> class Solver;

	//! Base Class for Evolutionary Algorithms
	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	class Solver_base
	{
	public:
		//! Solve wrapper function for Solvers, used for benchmarks
		std::vector<T> solver_bench(const std::string& problem_name);
	protected:
		//! Constructor
		Solver_base(const S<T>& i_solver_struct, const F& i_f, const C& i_c) :
			solver_struct{ i_solver_struct },
			individuals{ init_individuals() },
			min_cost{ individuals[0] },
			last_iter{ 0 },
			f{ i_f },
			c{ i_c },
			solved_flag{ false },
			timer{ 0 }
		{
			distribution = std::uniform_real_distribution<T>(0.0, 1.0);
			generator.discard(700000);
			for (auto& p : individuals)
			{
				p = randomise_individual();
				//! Check population constraints
				while (!c(p))
				{
					p = randomise_individual();
				}
			}
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
		size_t last_iter;
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
		//! Displays the results
		std::stringstream display_results();
		//! Write the results to a file
		void write_results_to_file(const std::string& problem_name);
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
		std::vector<std::vector<T>> individuals(solver_struct.npop, std::vector<T>(solver_struct.ndv));
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
	std::stringstream Solver_base<Derived, S, T, F, C>::display_results()
	{
		std::stringstream results;
		results << solver_struct.type << ",";
		if (!solved_flag)
		{
			results << "False" << ",";
		}
		else
		{
			results << "True" << ",";
		}
		results << min_cost << ",";
		results << f(min_cost) << ",";
		results << individuals.size() << ",";
		results << last_iter << ",";
		results << timer << ",";
		results << solver_struct.decision_variables << ",";
		results << solver_struct.stdev << ",";
		results << solver_struct.npop << ",";
		results << solver_struct.tol << ",";
		results << solver_struct.iter_max << ",";
		results << static_cast<Derived*>(this)->display_parameters().str() << "\n";
		return results;
	}

	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	void Solver_base<Derived, S, T, F, C>::write_results_to_file(const std::string& problem_name)
	{
		std::string filename;
		std::string stypes;
		stypes.append("Algorithm,Solved,Solution,Fitness,Population,Iterations,Elapsed Time,Initial Solution,Standard Deviation,Initial Population,Tolerance,Maximum Iterations");
		filename.append(problem_name);
		filename.append("-results.csv");
		std::ifstream in(filename);
		bool written_first_line = false;
		if (in.good())
		{
			std::string line;
			std::getline(in, line);
			if (line == stypes)
			{
				written_first_line = true;
			}
		}
		in.close();
		std::ofstream out;
		out.open(filename, std::ofstream::out | std::ofstream::app);
		if (!written_first_line)
		{
			out << stypes << "\n";
		}
		out << display_results().str();
	}

	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	std::vector<T> Solver_base<Derived, S, T, F, C>::solver_bench(const std::string& problem_name)
	{
		if (solver_struct.tol > std::abs(f(min_cost)))
		{
			timer = 0;
		}
		else
		{
			//! Time the computation
			const std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
			static_cast<Derived*>(this)->run_algo();
			const std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
			const std::chrono::duration<double> elapsed_seconds = end - start;
			timer = elapsed_seconds.count();
		}
		//! Return minimum cost individual
		if (solver_struct.print_to_output)
		{
			std::cout << display_results().str();
		}
		else {};
		if (solver_struct.print_to_file)
		{
			write_results_to_file(problem_name);
		}
		else {};
		return min_cost;
	}
	/*!
	/brief Solver wrapper function, interface to solvers : free function used for benchmarks
	*/
	template<typename F, typename C, template<typename> class S, typename T>
	std::vector<T> solve(const F& f, const C& c, const S<T>& solver_struct, const std::string& problem_name)
	{
		Solver<S, T, F, C> solver{ solver_struct, f, c };
		return solver.solver_bench(problem_name);
	}
}