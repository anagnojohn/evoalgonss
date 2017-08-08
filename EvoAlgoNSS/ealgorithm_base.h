#pragma once

#include <vector>
#include <assert.h>
#include <utility>
#include <chrono>
#include <ctime>
#include <random>
#include <vector>
#include <tuple>

template<typename T>
struct EAstruct
{
public:
	EAstruct(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: decision_variables{ i_decision_variables }, stdev{ i_stdev }, npop{ i_npop }, tol{ i_tol }, iter_max{ i_iter_max }, ndv { i_decision_variables.size() }
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
	// Decision Variables
	const std::vector<T> decision_variables;
	// Standard deviation of the decision variables
	const std::vector<T> stdev;
	// Size of the population
	const size_t npop;
	// Tolerance
	const T tol;
	// Number of maximum iterations
	const size_t iter_max;
	// Number of decision variables
	const size_t ndv;
};

template<typename T>
std::vector<T> randomise_individual(const std::vector<T>& decision_variables, const std::vector<T>& stdev)
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
std::vector<std::vector<T>> init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev)
{
	const auto& ndv = decision_variables.size();
	std::vector<std::vector<T>> individuals(npop, std::vector<T>(ndv));
	for (auto& p : individuals)
	{
		p = randomise_individual(decision_variables, stdev);
	}
	return individuals;
}

// Template Class for Solvers
template<typename T, typename S>
class Solver
{
public:
	template<typename F, typename C> Solver(const S& solver_struct, F f, C c)
	{

	}
};

// Base Class for Evolutionary Algorithms
template<typename T>
class Solver_base
{
public:
	template<typename F, typename C> Solver_base(const EAstruct<T>& solver_struct, F f, C c)
		: npop{ solver_struct.npop }, ndv{ solver_struct.ndv }, tol{ solver_struct.tol }, iter_max{ solver_struct.iter_max },
		individuals{ init_individuals(solver_struct.decision_variables, solver_struct.npop, solver_struct.stdev) }, iter{ 0 }
	{
		min_cost = individuals[0];
		decision_variables = solver_struct.decision_variables;
		stdev = solver_struct.stdev;
		population_constraints_checker(c);
		find_min_cost(f);
	}
	std::vector<T> ret_min_cost() { return min_cost; };
	T ret_min_cost() const { return min_cost; };
	size_t ret_npop() const { return npop; };
	T ret_fitness_cost() const { return fitness_cost; };
	size_t ret_iter() const { return iter; };
	T ret_tol() const { return tol; };
	template<typename F, typename C> void run_algo(F f, C c) {};
	const std::string type = "";
protected:
	// The decision variables
	std::vector<T> decision_variables;
	// The standard deviation of decision variables
	std::vector<T> stdev;
	// Size of the population
	size_t npop;
	// Tolerance
	T tol;
	// Number of maximum iterations
	size_t iter_max;
	// Number of decision variables
	size_t ndv;
	// Population
	std::vector<std::vector<T>> individuals;
	// Best fitness
	T fitness_cost;
	// Last iteration to solution
	size_t iter;
	// Best solution
	std::vector<T> min_cost;
	// Find best solution
	template<typename F> void find_min_cost(F f);
	// Check population constraints
	template<typename C> void population_constraints_checker(C c)
	{
		for (auto& p : individuals)
		{
			while (!c(p))
			{
				p = randomise_individual(decision_variables, stdev);
			}
		}
	};
};

// Find the minimum cost individual of the fitness function for the population
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

template<typename T, typename F, typename C, typename S>
std::vector<T> solve(F f, C c, const S& solver_struct)
{
	Solver<T, S> solver(solver_struct, f, c);
	std::cout << solver.type << " used as solver" << "\n";
	if (solver.ret_tol() > std::abs(solver.ret_fitness_cost()))
	{
		T timer = 0;
		std::cout << "Optimum solution: " << solver.ret_min_cost() << " Fitness Value: " << solver.ret_fitness_cost() << "\n";
		std::cout << "Population: " << solver.ret_npop() << " Solved at iteration: " << solver.ret_iter() << "\n";
		std::cout << "Elapsed time in seconds: " << timer << "\n";
	}
	// Time the computation
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	solver.run_algo(f, c);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	T timer = elapsed_seconds.count();
	std::cout << "Optimum solution: " << solver.ret_min_cost() << " Fitness Value: " << solver.ret_fitness_cost() << "\n";
	std::cout << "Population: " << solver.ret_npop() << " Solved at iteration: " << solver.ret_iter() << "\n";
	std::cout << "Elapsed time in seconds: " << timer << "\n";
	// Return minimum cost individual
	return solver.ret_min_cost();
}