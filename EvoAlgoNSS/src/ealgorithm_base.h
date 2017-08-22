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
#include "utilities.h"

//! Evolutionary Algorithms
namespace ea
{
	using namespace utilities;
	/** \struct EA_base
	*  \brief Evolutionary algorithm stucture base
	*/
	template<typename T>
	struct EA_base
	{
	public:
		/** \brief The floating point number type used for type deduction */
		using fp_type = T;
		/** \brief Initial Decision Variables */
		const std::vector<T> decision_variables;
		/** \brief Standard deviation of the decision variables */
		const std::vector<T> stdev;
		/** \brief Size of the population */
		const size_t npop;
		/** \brief Tolerance */
		const T tol;
		/** \brief Number of maximum iterations */
		const size_t iter_max;
		/** \brief Number of decision variables */
		const size_t ndv;
		/** \brief Use penalty function or not */
		const bool use_penalty_method;
		/** \brief Constraints type */
		const Constraints_type constraints_type;
		/** \brief Print to output or not */
		const bool print_to_output;
		/** \brief Print to file or not */
		const bool print_to_file;
	protected:
		/** \fn EA_base(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max,
			const bool& i_use_penalty_method, const Constraints_type& i_constraints_type, const bool& i_print_to_output, const bool& i_print_to_file)
		*	\brief Constructor
		*	\param i_decision_variables The starting values of the decision variables
		*	\param i_stdev The standard deviation
		*	\param i_npop The population size
		*	\param i_tol The tolerance
		*	\param i_iter_max The maximum number of iterations
		*	\param i_use_penalty_method Whether to used penalties or not
		*	\param i_constraints_type What kind of constraints to use
		*	\param i_print_to_output Whether to print to terminal or not
		*	\param i_print_to_file Whether to print to a file or not
		*	\return A EA_base<T> object
			*/
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
	/** \brief Random device/ Random number generator */
	std::random_device rd;
	/** \fn generator(rd())
	*  \brief Pseudo-random number generator
	*  \return A random number
	*/
	std::mt19937_64 generator(rd());
	/*! \class Solver
	*  \brief Template Class for Solvers */
	template<template<typename> class S, typename T, typename F, typename C> class Solver;
	
	/*! \class Solver_base
	*  \brief Base Class for Evolutionary Algorithms
	*  \details The fitness and constraints functions are copied
	*  so that even if they are not available in the current scope, the solver will still execute properly.
	*  At the same time, std::function could be have used, thus eliminating the need for template parameters F and C.
	*  However, that comes at a runtime cost, since calls to the functions would be virtual and there is a possibility
	*  that allocation could happen on the heap.*/
	template<typename Derived, template<typename> class S, typename T, typename F, typename C>
	class Solver_base
	{
	public:
		/*! \fn solver_bench(const std::string& problem_name)
		*  \brief Solve wrapper function for Solvers, used for benchmarks
		*  \param problem_name The name of the problem in std::string form
		*  \return The solution vector
		*/
		std::vector<T> solver_bench(const std::string& problem_name);
	protected:
		/*! \fn Solver_base(const S<T>& i_solver_struct, const F& i_f, const C& i_c)
		*  \brief Constructor
		*  \param i_solver_struct The parameter structure that is used to construct the solver
		*  \param i_f A reference to the objective function
		*  \param i_c A reference to the constraints function
		*  \return A Solver_base<Derived, S, T, F, C> object
		*/
		Solver_base(const S<T>& i_solver_struct, const F& i_f, const C& i_c) :
			solver_struct{ i_solver_struct },
			individuals{ init_individuals(i_solver_struct) },
			min_cost{ individuals[0] },
			last_iter{ 0 },
			f{ i_f },
			c{ i_c },
			solved_flag{ false },
			timer{ 0 },
			distribution{ std::uniform_real_distribution<T>(0.0, 1.0) }
		{
			generator.discard(700000);
			find_min_cost();
		}
		/** \brief Internal copy of the structure used for parameters of the algorithm */
		const S<T> solver_struct;
		/** \brief Population */
		std::vector<std::vector<T>> individuals;
		/** \brief Best solution / lowest fitness */
		std::vector<T> min_cost;
		/** \brief  Last iteration to solution */
		size_t last_iter;
		/** \brief Copy of the fitness function passed as a lambda */
		F f;
		/** \brief Copy of the constraints function passed as a lambda */
		C c;
        /** \brief A flag which determines if the solver has already solved the problem */
        bool solved_flag;
        /** \brief The timer used for benchmarks */
        T timer;
		/** \brief Uniform real distribution */
		std::uniform_real_distribution<T> distribution;
		/*! \fn randomise_individual()
		*  \brief Returns a randomised individual using the initial decision variables and standard deviation
		*  \return A randomised individual of type std::vector<T>, where T is a floating-point number type.
		*/
		std::vector<T> randomise_individual();
		/*! \fn init_individuals()
		*  \brief Initialises the population by randomising aroung the decision variables using the given standard deviation
		*  \return The population after checking the constraints of the optimisation problem
		*/
		std::vector<std::vector<T>> init_individuals(const S<T>& solver_struct);
		/*! \fn find_min_cost()
		*  \brief Find the minimum cost individual of the fitness function for the population
		*  \return void
		*/
		void find_min_cost();
		/*! \fn display_results()
		*  \brief Display the results of execution of an algorithm as well as its parameters
		*  \return A std::stringstream of the results
		*/
		std::stringstream display_results();
		/*! \fn write_results_to_file(const std::string& problem_name)
		*  \brief Write the results to a file
		*  \param problem_name The name of the problem in std::string form
		*  \return void
		*/
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
	std::vector<std::vector<T>> Solver_base<Derived, S, T, F, C>::init_individuals(const S<T>& solver_struct)
	{
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
	
	/*! \fn solve(const F& f, const C& c, const S<T>& solver_struct, const std::string& problem_name)
	*  \brief Solver wrapper function, interface to solvers : free function used for benchmarks
	*  \param f The objective function
	*  \param c The constraints function
	*  \param solver_struct The parameter structure of the solver
	*  \param problem_name The name of the problem in std::string form
	*  \return The solution vector
	*/
	template<typename F, typename C, template<typename> class S, typename T>
	std::vector<T> solve(const F& f, const C& c, const S<T>& solver_struct, const std::string& problem_name)
	{
		Solver<S, T, F, C> solver{ solver_struct, f, c };
		return solver.solver_bench(problem_name);
	}
}