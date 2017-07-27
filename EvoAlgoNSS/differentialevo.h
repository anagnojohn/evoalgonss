#pragma once

#include <random>
#include <vector>
#include <tuple>
#include <chrono>
#include <ctime>
#include "ealgorithm_base.h"

template<typename T>
struct DEstruct : EA_base<T>
{
public:
	DEstruct(const T& i_cr, const T& i_f_param, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: cr{ i_cr }, f_param{ i_f_param }, EA_base{ i_npop, i_tol, i_iter_max }
	{
		assert(cr > 0 && cr <= 1);
		assert(f_param > 0 && f_param <= 1);
	}
	// Crossover Rate
	T cr;
	// Mutation Scale Fuctor
	T f_param;
};

//Differential Evolution Algorithm Class
template<typename T, typename F>
class Solver<T, F, DEstruct<T>>
{
public:
	// Constructor for the Differential Evolution Class
	Solver(const DEstruct<T>& de, const Population<T>& popul)
	{
		cr = de.cr;
		f_param = de.f_param;
		npop = popul.npop;
		individuals = popul.individuals;
		ndv = popul.ndv;
		tol = de.tol;
		iter_max = de.iter_max;
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		for (auto i = 0; i < npop; ++i)
		{
			indices.push_back(i);
		}
		std::uniform_int_distribution<size_t> i_ind_distribution(0, indices.size() - 1);
		ind_distribution = i_ind_distribution;
	}
	std::tuple<std::vector<T>, T, size_t, double> solve(F f, const T& opt);
private:
	// Crossover Rate
	T cr;
	// Mutation Scale Fuctor
	T f_param;
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
	// Indices of population
	std::vector<size_t> indices;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	std::random_device generator;
	// Uniform Real Distribution used for generating the random vectors for mutation and epsilon for crossover 
	std::uniform_real_distribution<T> distribution;
	std::uniform_int_distribution<size_t> ind_distribution;
	void construct_donor(std::vector<T>& donor);
	void construct_trial(const std::vector<T>& target, const std::vector<T>& donor, std::vector<T>& trial);
};

// Method that constructs the donor vector
template<typename T, typename F>
void Solver<T, F, DEstruct<T>>::construct_donor(std::vector<T>& donor)
{
	std::vector<size_t> r_i;
	while (r_i.size() < 3)
	{
		r_i.push_back(indices[ind_distribution(engine)]);
		if (r_i.size() > 1 && r_i.end()[-1] == r_i.end()[-2])
		{
			r_i.pop_back();
		}
	}
	for (auto j = 0; j < ndv; ++j)
	{
		donor[j] = individuals[r_i[0]][j] + f_param * (individuals[r_i[1]][j] - individuals[r_i[2]][j]);
	}
}

// Method that constructs the trial vector
template<typename T, typename F>
void Solver<T, F, DEstruct<T>>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor, std::vector<T>& trial)
{
	for (auto j = 0; j < ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = indices[ind_distribution(engine)];
		if (epsilon <= cr || j == jrand)
		{
			trial[j] = donor[j];
		}
		else
		{
			trial[j] = target[j];
		}
	}
}

// Solve function overload for the Differential Evolution class
template<typename T, typename F>
std::tuple<std::vector<T>, T, size_t, double> Solver<T,F, DEstruct<T>>::solve(F f, const T& opt)
{
	std::cout << "Differential Evolution used as solver" << "\n";
	// Find the minimum cost individual of the fitness function for the population
	std::vector<T> min_cost = individuals[0];
	std::vector<T> donor = min_cost;
	std::vector<T> trial = donor;
	for (const auto& p : individuals)
	{
		if (f(min_cost) > f(p))
		{
			min_cost = p;
		}
	}
	T fitness_cost = f(min_cost);
	size_t last_iter = 0;
	// Time the computation
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	// Differential Evolution starts here
	for (auto iter = 0; iter < iter_max; ++iter)
	{
		// Stopping Criteria
		if (tol > std::abs(fitness_cost - opt))
		{
			last_iter = iter;
			break;
		}
		for (auto i = 0; i < npop; ++i)
		{
			// Construct donor and trial vectors
			construct_donor(donor);
			construct_trial(individuals[i], donor, trial);
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
		fitness_cost = f(min_cost);
		last_iter = iter;
	}
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	T timer = elapsed_seconds.count();
	// Return minimum cost individual
	return { min_cost, fitness_cost, last_iter, timer };
}