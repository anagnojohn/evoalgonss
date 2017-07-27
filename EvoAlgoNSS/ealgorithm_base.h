#pragma once

#include <vector>
#include <assert.h>

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
	Solver(const S& solver_struct, const Population<T>& popul)
	{

	}
};

// Base Class for Evolutionary Algorithms
template<typename T, typename F>
class Solver<T,F, EAstruct<T>>
{
public:
	Solver(const EAstruct<T>& solver_struct, const Population<T>& popul)
	{
		npop = popul.npop;
		individuals = popul.individuals;
		ndv = popul.ndv;
		tol = solver_struct.tol;
		iter_max = solver_struct.iter_max;
	}
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
};

template<typename T, typename F, typename S>
Solver<T, F, S> create_solver(const S& solver_struct, const Population<T>& popul)
{
	Solver<T, F, S> solver(solver_struct, popul);
	return solver;
}