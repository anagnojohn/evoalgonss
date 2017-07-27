#pragma once

#include <random>
#include <vector>
#include <tuple>
#include <chrono>
#include <ctime>
#include "ealgorithm_base.h"

template<typename T>
T euclid_distance(const std::vector<T>& x, const std::vector<T>& y)
{
	T sum = 0;
	for (auto i = 0; i < x.size(); ++i)
	{
		sum = sum + std::pow(x[i] - y[i], 2);
	}
	return std::sqrt(sum);
}

template<typename T>
struct PSOstruct : EA_base<T>
{
public:
	PSOstruct(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: c1{ i_c1 }, c2{ i_c2 }, sneigh{ i_sneigh }, EA_base{ i_npop, i_tol, i_iter_max }
	{
		assert(c1 > 0);
		assert(c2 > 0);
		assert(sneigh > 0);
		assert(sneigh < npop);
	}
	T c1;
	T c2;
	// Size of each neighbourhood
	size_t sneigh;
};

// Local Best Particle Swarm Optimisation Class
template<typename T, typename F>
class Solver<T, F, PSOstruct<T>>
{
public:
	Solver(const PSOstruct<T>& pso, const Population<T>& popul)
	{
		c1 = pso.c1;
		c2 = pso.c2;
		sneigh = pso.sneigh;
		npop = popul.npop;
		individuals = popul.individuals;
		ndv = popul.ndv;
		tol = pso.tol;
		iter_max = pso.iter_max;
		init_pso();
	}
	std::tuple<std::vector<T>, T, size_t, double> solve(F f, const T& opt);
protected:
	T c1;
	T c2;
	// Neighbourhood size
	size_t sneigh;
	// Size of the population
	size_t npop;
	// Tolerance
	T tol;
	// Number of maximum iterations
	size_t iter_max;
	// Population
	std::vector<std::vector<T>> individuals;
	std::vector<std::vector<T>> personal_best;
	std::vector<std::vector<T>> local_best;
	std::vector<std::vector<T>> velocity;
	// Number of neighbourhoods
	size_t nneigh;
	// Neighbourhoods
	std::vector<std::vector<size_t>> neighbourhoods;
	// Indices of the population
	std::vector<size_t> indices;
	// Number of decision variables
	size_t ndv;
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
	void init_pso();
	void velocity_update(const size_t& iter);
	void set_neighbourhoods();
	void set_local_best(F f);
	std::vector<std::vector<T>> generate_r();
	std::vector<T> find_min_local_best(F f);
};

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::init_pso()
{
	std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
	distribution = i_distribution;
	nneigh = static_cast<size_t>(std::ceil(npop / sneigh));
	neighbourhoods.resize(nneigh);
	velocity.resize(npop, std::vector<T>(ndv));
	local_best.resize(nneigh);
	for (auto& p : velocity)
	{
		for (auto& n : p)
		{
			n = 0.0;
		}
	}
	for (auto k = 0; k < neighbourhoods.size(); ++k)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			local_best[k].push_back(individuals[0][j]);
		}
	}
	for (const auto& p : individuals)
	{
		personal_best.push_back(p);
	}
	for (auto i = 0; i < individuals.size(); ++i)
	{
		indices.push_back(i);
	}
	set_neighbourhoods();
}

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::set_neighbourhoods()
{
	const auto& nneigh = neighbourhoods.size();
	size_t counter = 0;
	for (auto k = 0; k < nneigh; ++k)
	{
		while (counter < sneigh)
		{
			neighbourhoods[k].push_back(indices[k * sneigh + counter]);
			counter = counter + 1;
		}
	}
}

template<typename T, typename F>
std::vector<std::vector<T>> Solver<T, F, PSOstruct<T>>::generate_r()
{
	std::vector<std::vector<T>> r(2, std::vector<T>(ndv));
	for (auto i = 0; i < 2; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			r[i][j] = (distribution(generator));
		}
	}
	return r;
}

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::velocity_update(const size_t& iter)
{
	const auto& r = generate_r();
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			for (auto j = 0; j < ndv; ++j)
			{
				{
					velocity[i][j] = velocity[i][j] + c1 * r[0][j] * (personal_best[i][j]
						- individuals[i][j]) + c2 * r[1][j] * (local_best[k][j] - individuals[i][j]);
					individuals[i][j] = individuals[i][j] + velocity[i][j];
				}
			}
		}
	}
}

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::set_local_best(F f)
{
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			if (f(personal_best[i]) < f(local_best[k]))
			{
				for (auto j = 0; j < ndv; j++)
				{
					local_best[k][j] = personal_best[i][j];
				}
			}
		}
	}
}


template<typename T, typename F>
std::vector<T> Solver<T, F, PSOstruct<T>>::find_min_local_best(F f)
{
	std::vector<T> min_cost(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		min_cost[j] = local_best[0][j];
	}
	for (auto k = 0; k < nneigh; ++k)
	{
		if (f(local_best[k]) < f(min_cost))
		{
			min_cost = local_best[k];
		}
	}
	return min_cost;
}

template<typename T, typename F>
std::tuple<std::vector<T>, T, size_t, double> Solver<T, F, PSOstruct<T>>::solve(F f, const T& opt)
{
	std::cout << "Local Best Particle Swarm used as solver" << "\n";
	T rmax;
	std::vector<T> distance(npop);
	std::vector<T> min_cost(ndv);
	set_local_best(f);
	min_cost = find_min_local_best(f);
	T fitness_cost = f(min_cost);
	size_t last_iter = 0;
	// Time the computation
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	// Local Best Particle Swarm starts here
	for (auto iter = 0; iter < iter_max; ++iter)
	{	
		for (auto i = 0; i < npop; ++i)
		{
			distance[i] = euclid_distance(individuals[i], min_cost);
		}
		rmax = distance[0];
		for (auto i = 0; i < npop; ++i)
		{
			if (rmax < distance[i])
			{
				rmax = distance[i];
			}
		}
		if (tol > std::abs(fitness_cost - opt) || rmax < tol)
		{
			last_iter = iter;
			break;
		}
		velocity_update(iter);
		for (auto i = 0; i < npop; ++i)
		{
			if (f(individuals[i]) < f(personal_best[i]))
			{
				personal_best[i] = individuals[i];
			}
		}
		for (auto k = 0; k < neighbourhoods.size(); ++k)
		{
			for (auto l = 0; l < neighbourhoods[k].size(); ++l)
			{
				auto i = neighbourhoods[k][l];
				if (f(personal_best[i]) < f(local_best[k]))
				{
						local_best[k] = personal_best[i];
				}
			}
		}
		min_cost = find_min_local_best(f);
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