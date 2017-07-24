#pragma once

#include <random>
#include <vector>
#include "ealgorithm_base.h"

//Differential Evolution Class
template<typename T>
class DifferentialEvo
{
public:
	// Constructor for the Differential Evolution Class
	DifferentialEvo(const T& i_cr, const T& i_f_param)
		: cr(i_cr), f_param(i_f_param)
	{
		assert((cr > 0 && cr <= 1) && (f_param > 0 && f_param <= 1));
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	std::vector<T> construct_donor(const std::vector< std::vector<T> >& individuals);
	std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
	// Setters
	void set_cr(const T& cr) { assert(cr > 0 && cr <= 1); this->cr = cr; };
	void set_f_param(const T& f_param) { assert(f_param > 0 && f_param <= 1); this->f_param = f_param; };
private:
	// Crossover Rate
	const T cr;
	// Mutation Scale Fuctor
	const T f_param;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	std::random_device generator;
	// Uniform Real Distribution used for generating the random vectors for mutation and epsilon for crossover 
	std::uniform_real_distribution<T> distribution;
};

// Method that constructs the donor vector
template<typename T>
std::vector<T> DifferentialEvo<T>::construct_donor(const std::vector< std::vector<T> >& individuals)
{
	const auto& ndv = individuals[0].size();
	const auto& npop = individuals.size();
	std::vector<T> donor(ndv);
	std::vector<size_t> r_i;
	std::vector<size_t> indices;
	for (auto i = 0; i < npop; ++i)
	{
		indices.push_back(i);
	}
	std::uniform_int_distribution<size_t> distribution(0, indices.size() - 1);
	while (r_i.size() < 4)
	{
		r_i.push_back(indices[distribution(engine)]);
		if (r_i.size() > 1 && r_i.end()[-1] == r_i.end()[-2])
		{
			r_i.pop_back();
		}
	}
	for (auto j = 0; j < ndv; ++j)
	{
		donor[j] = individuals[r_i[0]][j] + f_param * (individuals[r_i[1]][j] - individuals[r_i[2]][j]);
	}
	return donor;
}

// Method that constructs the trial vector
template<typename T>
std::vector<T> DifferentialEvo<T>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
{
	auto ndv = target.size();
	std::vector<T> trial(ndv);
	std::vector<size_t> indices;
	for (auto i = 0; i < ndv; ++i)
	{
		indices.push_back(i);
	}
	std::uniform_int_distribution<size_t> size_distribution(0, indices.size() - 1);
	for (auto j = 0; j < ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = indices[size_distribution(engine)];
		if (epsilon <= cr || j == jrand)
		{
			trial[j] = donor[j];
		}
		else
		{
			trial[j] = target[j];
		}
	}
	return trial;
}

// Solve function overload for the Differential Evolution class
template<typename T, typename F>
std::vector<T> solve(F f, const T& opt, DifferentialEvo<T>& de, const EAparams<T>& ea)
{
	// Find the minimum cost individual of the fitness function for the population
	const auto& tol = ea.get_tol();
	const auto& iter_max = ea.get_iter_max();
	const auto& npop = ea.get_npop();
	const auto& ndv = ea.get_ndv();
	const auto& stdev = ea.get_stdev();
	auto individuals = ea.get_individuals();
	std::vector<T> min_cost = individuals[0];
	for (const auto& p : individuals)
	{
		if (f(min_cost) > f(p))
		{
			min_cost = p;
		}
	}
	// Differential Evolution starts here
	for (auto g = 0; g < iter_max; ++g)
	{
		// Stopping Criteria
		if (tol > std::abs(f(min_cost) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		for (auto i = 0; i < npop; ++i)
		{
			// Construct donor and trial vectors
			std::vector<T> donor = de.construct_donor(individuals);
			std::vector<T> trial = de.construct_trial(individuals[i], donor);
			// Replace individual i with trial if trial's cost is lower / Selection step
			if (f(trial) <= f(individuals[i]))
			{
				individuals[i] = trial;
			}
		}
		// Recalculate minimum cost individual of the population
		for (const auto& p : individuals)
		{
			if (f(min_cost) > f(p))
			{
				min_cost = p;
			}
		}
	}
	T error = std::abs(f(min_cost) - opt);
	// Return minimum cost individual
	return min_cost;
}