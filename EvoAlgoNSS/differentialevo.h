#pragma once

#include "ealgorithm_base.h"

//Differential Evolution Class
template<typename T>
class DifferentialEvo : public EA_base<T>
{
public:
	// Constructor for the Differential Evolution Class
	DifferentialEvo(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max, const T& i_cr, const T& i_f_param, const std::vector<T>& stdev)
		: cr(i_cr), f_param(i_f_param)
	{
		set_solver(decision_variables, npop, tol, iter_max);
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	std::vector<T> construct_donor();
	std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
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
std::vector<T> DifferentialEvo<T>::construct_donor()
{
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
std::vector<T> solve(F f, const T& opt, DifferentialEvo<T>& de)
{
	// Find the minimum cost individual of the fitness function for the population
	const auto& tol = de.get_tol();
	const auto& iter_max = de.get_iter_max();
	const auto& npop = de.get_npop();
	const auto& ndv = de.get_ndv();
	const auto& stdev = de.get_stdev();
	const auto& individuals = de.get_individuals();
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
		for (auto i = 0; i < individuals.size(); ++i)
		{
			// Construct donor and trial vectors
			std::vector<T> donor = de.construct_donor();
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