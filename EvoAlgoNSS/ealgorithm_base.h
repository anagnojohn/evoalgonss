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
struct Population
{
public:
	Population(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop)
		: decision_variables{ i_decision_variables }, stdev{ i_stdev }, npop{ i_npop }, ndv{ i_decision_variables.size() }, 
		individuals { init_individuals(i_decision_variables, i_npop, i_stdev) }
	{
		assert(decision_variables.size() > 0);
		assert(decision_variables.size() == stdev.size());
		for (const auto& p : stdev)
		{
			assert(p > 0);
		}
		assert(npop > 0);
	}
	// Decision Variables
	const std::vector<T> decision_variables;
	// Standard deviation of the decision variables
	const std::vector<T> stdev;
	// Size of the population
	const size_t npop;
	// Population
	const std::vector<std::vector<T>> individuals;
	// Number of decision variables
	const size_t ndv;
};

template<typename T>
struct EAstruct
{
public:
	EAstruct(const size_t& i_npop, const T& i_tol, const size_t& i_iter_max) : npop{ i_npop }, tol{ i_tol }, iter_max{ i_iter_max }
	{
		assert(npop > 0);
		assert(tol > 0);
		assert(iter_max > 0);
	}
	// Size of the population
	const size_t npop;
	// Tolerance
	const T tol;
	// Number of maximum iterations
	const size_t iter_max;
};

template<typename T>
std::vector<std::vector<T>> init_individuals(const std::vector<T>& decision_variables, const size_t& npop, const std::vector<T>& stdev)
{
	std::random_device generator;
	T epsilon;
	const auto& ndv = decision_variables.size();
	std::vector<std::vector<T>> individuals;
	individuals.resize(npop, std::vector<T>(ndv));
	for (auto& p : individuals)
	{
		p = decision_variables;
	}
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			std::normal_distribution<T> ndistribution(0, stdev[j]);
			epsilon = ndistribution(generator);
			individuals[i][j] = individuals[i][j] + epsilon;
		}
	}
	return individuals;
}

// Template Class for Solvers
template<typename T, typename F, typename S>
class Solver
{
public:
	Solver(const S& solver_struct, const Population<T>& popul, F f)
	{

	}
};

// Base Class for Evolutionary Algorithms
template<typename T, typename F>
class Solver<T, F, EAstruct<T>>
{
public:
	Solver(const EAstruct<T>& solver_struct, const Population<T>& popul) : npop{ popul.npop },
		individuals{ popul.individuals }, ndv{ popul.ndv }, tol{ solver_struct.tol }, iter_max{ solver_struct.iter_max }
	{
		min_cost = individuals[0];
	}
	std::tuple<std::vector<T>, T, size_t, double> solve(F f, const T& opt);
protected:
	// Size of the population
	size_t npop;
	// Tolerance
	T tol;
	// Number of maximum iterations
	size_t iter_max;
	// Population
	std::vector<std::vector<T>> individuals;
	// Number of decision variables
	size_t ndv;
	// Best fitness
	T fitness_cost;
	// Last iteration to solution
	size_t last_iter;
	// Best solution
	std::vector<T> min_cost;
	// Find best solution
	void find_min_cost(F f);
	// Get the type of solver
	virtual std::string get_type_of_solver() { return ""; };
	// Run a solver
	virtual void run_algo(F f, const T& opt) {};
};

template<typename T, typename F>
void Solver<T, F, EAstruct<T>>::find_min_cost(F f)
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

template<typename T, typename F>
std::tuple<std::vector<T>, T, size_t, double> Solver<T, F, EAstruct<T>>::solve(F f, const T& opt)
{
	std::cout << get_type_of_solver() << " used as solver" << "\n";
	// Find the minimum cost individual of the fitness function for the population
	find_min_cost(f);
	// Time the computation
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	run_algo(f, opt);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	T timer = elapsed_seconds.count();
	// Return minimum cost individual
	return { min_cost, fitness_cost, last_iter, timer };
}

template<typename T, typename F, typename S, typename C>
Solver<T, F, S> create_solver(const S& solver_struct, const Population<T>& popul)
{
	Solver<T, F, S> solver(solver_struct, popul);
	return solver;
}