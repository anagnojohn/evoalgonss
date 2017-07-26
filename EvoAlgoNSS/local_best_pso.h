#pragma once

#include <random>
#include <vector>
#include <tuple>
#include <chrono>
#include <ctime>
#include "ealgorithm_base.h"

template<typename T>
class LocalBestPSO : public EA_base<T>
{
public:
	LocalBestPSO(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: c1{ i_c1 }, c2{ i_c2 }, sneigh{ i_sneigh }, EA_base{ i_npop, i_tol, i_iter_max }
		
	{
		assert(c1 > 0);
		assert(c2 > 0);
		assert(sneigh > 0);
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		vmax_flag = false;
		intertia_flag = false;
	}
	LocalBestPSO(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const T& i_w, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: c1{ i_c1 }, c2{ i_c2 }, sneigh{ i_sneigh }, w{ i_w }, EA_base{ i_npop, i_tol, i_iter_max }
	{
		assert(c1 > 0);
		assert(c2 > 0);
		assert(sneigh > 0);
		assert(w > 0);
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		vmax_flag = false;
		inertia_flag = true;
	}
	LocalBestPSO(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const std::vector<T>& i_vmax, const T& i_alpha, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: c1{ i_c1 }, c2{ i_c2 }, sneigh{ i_sneigh }, alpha{ i_alpha }, vmax{ i_vmax }, EA_base{ i_npop, i_tol, i_iter_max }
	{
		assert(c1 > 0);
		assert(c2 > 0);
		assert(sneigh > 0);
		assert(alpha > 0);
		for (const auto& p : vmax) { assert(p > 0); };
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		vmax_flag = true;
		inertia_flag = false;
	}
	void init_pso(const std::vector<std::vector<T>>& individuals, std::vector<std::vector<T>>& personal_best,
		std::vector<std::vector<T>>& local_best, std::vector<std::vector<T>>& velocity, std::vector<std::vector<size_t>>& neighbourhoods);
	void velocity_update(std::vector< std::vector<T> >& individuals, std::vector<std::vector<T>>& personal_best,
		std::vector<std::vector<T>>& local_best, std::vector< std::vector<T> >& velocity, const size_t& iter, const std::vector<std::vector<size_t>>& neighbourhoods);
	// Getters
	size_t get_neighbourhood_size() const { return sneigh; };
	// Setters
	void set_c1(const T& c1) { assert(c1 > 0); this->c1 = c1; };
	void set_c2(const T& c2) { assert(c2 > 0); this->c2 = c2; };
	void set_alpha(const T& alpha) { assert(alpha > 0); this->alpha = alpha; };
	void set_vmax(const std::vector<T>& vmax)
	{
		for (const auto& p : vmax)
		{ 
			assert(p > 0);
		}
		vmax_flag = true;
		inertia_flag = false;
		this->vmax = vmax;
	};
	void set_w(const T& w) { assert(w > 0); vmax_flag = false; inertia_flag = true; this->w = w; };
	void set_sneigh(const size_t& sneigh) { assert(sneigh > 0); this->sneigh = sneigh; };
private:
	const T c1;
	const T c2;
	const T alpha;
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
	T w;
	std::vector<T> vmax;
	bool vmax_flag;
	bool inertia_flag;
	void set_neighbourhoods(std::vector<std::vector<size_t>>& neighbourhoods, std::vector<size_t>& indices);
	// Neighbourhood size
	const size_t sneigh;
	std::vector<std::vector<T>> generate_r(const size_t& ndv);
};

template<typename T>
void LocalBestPSO<T>::init_pso(const std::vector<std::vector<T>>& individuals, std::vector<std::vector<T>>& personal_best,
std::vector<std::vector<T>>& local_best, std::vector<std::vector<T>>& velocity, std::vector<std::vector<size_t>>& neighbourhoods)
{
	const auto& ndv = individuals[0].size();
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
	std::vector<size_t> indices;
	for (auto i = 0; i < individuals.size(); ++i)
	{
		indices.push_back(i);
	}
	set_neighbourhoods(neighbourhoods, indices);
}

template<typename T>
void LocalBestPSO<T>::set_neighbourhoods(std::vector<std::vector<size_t>>& neighbourhoods, std::vector<size_t>& indices)
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



template<typename T>
std::vector<std::vector<T>>LocalBestPSO<T>::generate_r(const size_t& ndv)
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

template<typename T>
void LocalBestPSO<T>::velocity_update(std::vector< std::vector<T> >& individuals, std::vector<std::vector<T>>& personal_best,
	std::vector<std::vector<T>>& local_best, std::vector< std::vector<T> >& velocity, const size_t& iter, const std::vector<std::vector<size_t>>& neighbourhoods)
{
	const auto& nneigh = neighbourhoods.size();
	const auto& ndv = individuals[0].size();
	const auto& r = generate_r(ndv);
	if (inertia_flag == false && vmax_flag == false)
	{
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
	if (vmax_flag)
	{
		assert(vmax.size() == ndv);
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
						if (velocity[i][j] > vmax[j])
						{
							velocity[i][j] = vmax[j];
						}
						individuals[i][j] = individuals[i][j] + velocity[i][j];
					}
				}
			}
		}
		for (auto& p : vmax)
		{
			p = (1 - std::pow(iter / iter_max, alpha)) * p;
		}
	}
	if (inertia_flag)
	{
		for (auto k = 0; k < nneigh; ++k)
		{
			for (auto l = 0; l < neighbourhoods[k].size(); ++l)
			{
				auto i = neighbourhoods[k][l];
				for (auto j = 0; j < ndv; ++j)
				{
					{
						velocity[i][j] = w * velocity[i][j] + c1 * r[0][j] * (personal_best[i][j]
							- individuals[i][j]) + c2 * r[1][j] * (local_best[k][j] - individuals[i][j]);
						individuals[i][j] = individuals[i][j] + velocity[i][j];
					}
				}
			}
		}
		w = ((w - 0.4) * (iter_max - iter)) / (iter_max + 0.4);
	}
}

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

template<typename T, typename F>
void set_local_best(F f, std::vector<std::vector<T>>& local_best, const std::vector<std::vector<size_t>>& neighbourhoods, const std::vector<std::vector<T>>& personal_best)
{
	const auto& nneigh = neighbourhoods.size();
	const auto& ndv = personal_best[0].size();
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
std::vector<T> find_min_local_best(const std::vector<std::vector<T>>& local_best, F f)
{
	const auto& nneigh = local_best.size();
	const auto& ndv = local_best[0].size();
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
std::tuple<std::vector<T>, T, size_t, double> solve(F f, const T& opt, const LocalBestPSO<T>& pso, const EAparams<T>& ea)
{
	std::cout << "Local Best Particle Swarm used as solver" << "\n";
	const auto& tol = pso.get_tol();
	const auto& iter_max = pso.get_iter_max();
	const auto& npop = ea.get_npop();
	const auto& ndv = ea.get_ndv();
	const auto& stdev = ea.get_stdev();
	auto individuals = ea.get_individuals();
	// Size of each neighbourhood
	const auto& sneigh = pso.get_neighbourhood_size();
	assert(sneigh <= individuals.size());
	// Number of neighbourhoods
	size_t nneigh = static_cast<size_t>(std::ceil(npop / sneigh));
	// Neighbourhoods
	std::vector<std::vector<size_t>> neighbourhoods(nneigh);
	std::vector<std::vector<T>> velocity(npop, std::vector<T>(ndv));
	std::vector<std::vector<T>> personal_best;
	std::vector<std::vector<T>> local_best(nneigh);
	pso.init_pso(individuals, personal_best, local_best, velocity, neighbourhoods);
	T rmax;
	std::vector<T> distance(npop);
	std::vector<T> min_cost(ndv);
	min_cost = find_min_local_best(local_best, f);
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
		pso.velocity_update(individuals, personal_best, local_best, velocity, iter, neighbourhoods);
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
		min_cost = find_min_local_best(local_best, f);
		fitness_cost = = f(min_cost);
		last_iter = iter;
	}
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	T timer = elapsed_seconds.count();
	// Return minimum cost individual
	return { min_cost, fitness_cost, last_iter, timer };
}