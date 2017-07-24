#pragma once

#include "differentialevo.h"

template<typename T, typename F, typename S = DifferentialEvo<T>>
class Solver
{
public:
	// Constructor for the Differential Evolution Class
	Solver(const DifferentialEvo<T>& i_de)
		: de(i_de)
	{
		set_solver();
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	std::vector<T> solve(F f, const T& opt);
private:
	void set_solver()
	{
		if (individuals.size() > 0)
		{
			individuals.erase(individuals.begin(), individuals.end());
		}
		create_individuals(de.decision_variables);
		init_epsilon();
	}
	void init_epsilon()
	{
		std::random_device generator;
		std::vector<T> diff(de.npop);
		//std::vector<T> sum(ndv);
		//std::vector<T> mean(ndv);
		//std::vector<T> stdev(ndv);
		//std::vector<T> epsilon(ndv);
		T epsilon;
		/*
		for (auto& p : sum)
		{
		p = 0.0;
		}
		for (auto& p : individuals)
		{
		for (auto j = 0; j < ndv; ++j)
		{
		sum[j] = sum[j] + p[j];
		}
		}
		for (auto j = 0; j < ndv; ++j)
		{
		mean[j] = sum[j] / (npop);
		}
		for (auto j = 0; j < ndv; ++j)
		{
		for (auto i = 0; i < npop; ++i)
		{
		diff[i] = individuals[i][j] - mean[j];
		}
		T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		stdev[j] = std::sqrt(sq_sum / npop);
		}
		*/
		for (auto i = 0; i < de.npop; ++i)
		{
			for (auto j = 0; j < de.ndv; ++j)
			{
				std::normal_distribution<> ndistribution(0, de.stdev[j]);
				epsilon = ndistribution(generator);
				individuals[i][j] = individuals[i][j] + epsilon;
			}
		}
	}

	void create_individuals(const std::vector<T>& decision_variables)
	{
		individuals.resize(de.npop, std::vector<T>(de.ndv));
		for (auto i = 0; i < de.npop; ++i)
		{
			individuals[i] = decision_variables;
		}
	}
	std::vector<T> construct_donor();
	std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	std::random_device generator;
	std::vector<std::vector<T>> individuals;
	// Uniform Real Distribution used for generating the random vectors for mutation and epsilon for crossover 
	std::uniform_real_distribution<T> distribution;
	const DifferentialEvo<T> de;
};

// Method that constructs the donor vector
template<typename T, typename F, typename S = DifferentialEvo<T>>
std::vector<T> Solver<T, F, S>::construct_donor()
{
	std::vector<T> donor(de.ndv);
	std::vector<size_t> r_i;
	std::vector<size_t> indices;
	for (auto i = 0; i < de.npop; ++i)
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
	for (auto j = 0; j < de.ndv; ++j)
	{
		donor[j] = individuals[r_i[0]][j] + de.f_param * (individuals[r_i[1]][j] - individuals[r_i[2]][j]);
	}
	return donor;
}

// Method that constructs the trial vector
template<typename T, typename F, typename S = DifferentialEvo<T>>
std::vector<T> Solver<T, F, S>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
{
	std::vector<T> trial(de.ndv);
	std::vector<size_t> indices;
	for (auto i = 0; i < de.ndv; ++i)
	{
		indices.push_back(i);
	}
	std::uniform_int_distribution<size_t> size_distribution(0, indices.size() - 1);
	for (auto j = 0; j < de.ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = indices[size_distribution(engine)];
		if (epsilon <= de.cr || j == jrand)
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
template<typename T, typename F, typename S = DifferentialEvo<T>>
std::vector<T> Solver<T, F, S>::solve(F f, const T& opt)
{
	// Find the minimum cost individual of the fitness function for the population
	std::vector<T> min_cost = individuals[0];
	for (const auto& p : individuals)
	{
		if (f(min_cost) > f(p))
		{
			min_cost = p;
		}
	}
	// Differential Evolution starts here
	for (auto g = 0; g < de.iter_max; ++g)
	{
		// Stopping Criteria
		if (de.tol > std::abs(f(min_cost) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		for (auto i = 0; i < individuals.size(); ++i)
		{
			// Construct donor and trial vectors
			std::vector<T> donor = construct_donor();
			std::vector<T> trial = construct_trial(individuals[i], donor);
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