#pragma once

#include "dependencies.h"

template<typename T>
class LocalBestPSO
{
public:
	LocalBestPSO(const std::vector<T>& i_decision_variables, const size_t& i_npop, const T& i_tol,const size_t& i_iter_max,
		const T& i_c1, const T& i_c2, const T& i_alpha, const T& i_w, const std::vector<T>& i_vmax,
		const bool& i_vmax_flag, const bool& i_inertia_flag, const std::vector<T>& i_stdev)
		: c1(i_c1), c2(i_c2), alpha(i_alpha), w(i_w), vmax(i_vmax), vmax_flag(i_vmax_flag), 
		inertia_flag(i_inertia_flag), tol(i_tol), iter_max(i_iter_max), stdev(i_stdev), npop(i_npop), ndv(i_decision_variables.size())
	{
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		particles = create_individuals(npop, i_decision_variables);
		init_epsilon(particles, stdev);
		nneigh = static_cast<size_t>(std::ceil(npop / 5));
		neighbourhoods.resize(nneigh);
		velocity.resize(npop, std::vector<T>(ndv));
		for (auto& p : velocity)
		{
			for (auto& n : p)
			{
				n = 0.0;
			}
		}
		for (auto k = 0; k < nneigh; ++k)
		{
			for (auto j = 0; j < ndv; ++j)
			{
				local_best[k].push_back(particles[0][j]);
			}
		}
		for (const auto& p : particles)
		{
			personal_best.push_back(p);
		}
		if (vmax_flag && inertia_flag)
		{
			vmax_flag = false;
		}
	}
	T euclid_distance(const std::vector<T>& x, const std::vector<T>& y);
	std::vector<std::vector<T>> generate_r();
	void velocity_update(std::vector < std::vector<T> >& r, size_t iter);
	void create_neighbourhoods();
	std::vector<std::vector<T>> particles;
	const size_t iter_max;
	const T tol;
	const std::vector<T> stdev;
	size_t npop;
	const size_t ndv;
	size_t nneigh;
	std::vector<std::vector<size_t>> neighbourhoods;
	std::vector<std::vector<T>> velocity;
	std::vector<std::vector<T>> personal_best;
	std::vector<std::vector<T>> local_best;
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
};

template<typename T>
std::vector<std::vector<T>>LocalBestPSO<T>::generate_r()
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
T LocalBestPSO<T>::euclid_distance(const std::vector<T>& x, const std::vector<T>& y)
{
	T sum = 0;
	for (auto i = 0; i < x.size(); ++i)
	{
		sum = sum + std::pow(x[i] - y[i], 2);
	}
	return std::sqrt(sum);
}

template<typename T>
void LocalBestPSO<T>::create_neighbourhoods()
{
	std::vector<size_t> indices;
	for (auto i = 0; i < npop; ++i)
	{
		indices.push_back(i);
	}
	size_t counter = 0;
	for (auto k = 0; k < nneigh; ++k)
	{
		while (counter < 5)
		{
			neighbourhoods[k].push_back(indices[k * 5 + counter]);
			counter = counter + 1;
		}
	}
}

template<typename T>
void LocalBestPSO<T>::velocity_update(std::vector < std::vector<T> >& r, size_t iter)
{
	if (!inertia_flag && !vmax_flag)
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
							- particles[i][j]) + c2 * r[1][j] * (local_best[k][j] - particles[i][j]);
					}
				}
			}
		}
	}
	if (vmax_flag)
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
							- particles[i][j]) + c2 * r[1][j] * (local_best[k][j] - particles[i][j]);
						if (velocity[i][j] > vmax[j])
						{
							velocity[i][j] = vmax[j];
						}
						particles[i][j] = particles[i][j] + velocity[i][j];
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
							- particles[i][j]) + c2 * r[1][j] * (local_best[k][j] - particles[i][j]);
					}
				}
			}
		}
		w = ((w - 0.4) * (iter_max - iter)) / (iter_max + 0.4);
	}
}


template<typename T, typename F>
std::vector<T> find_min_local_best(const std::vector<std::vector<T>>& local_best, F f)
{
	auto nneigh = local_best.size();
	auto ndv = local_best[0].size();
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
std::vector<T> solve(F f, const T& opt, LocalBestPSO<T>& pso)
{
	T rmax;
	std::vector<T> distance;
	std::vector<T> min_cost(pso.ndv);
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	pso.create_neighbourhoods();
	for (auto k = 0; k < pso.nneigh; ++k)
	{
		for (auto l = 0; l < pso.neighbourhoods[k].size(); ++l)
		{
			auto i = pso.neighbourhoods[k][l];
			if (f(pso.personal_best[i]) < f(pso.local_best[k]))
			{
				for (auto j = 0; j < pso.ndv; j++)
				{
					pso.local_best[k][j] = pso.personal_best[i][j];
				}
			}
		}
	}
	min_cost = find_min_local_best(pso.local_best, f);
	for (auto iter = 0; iter < pso.iter_max; ++iter)
	{	
		for (auto i = 0; i < pso.npop; ++i)
		{
			distance.push_back(pso.euclid_distance(pso.particles[i], min_cost));
		}
		rmax = distance[0];
		for (auto i = 0; i < pso.npop; ++i)
		{
			if (rmax < distance[i])
			{
				rmax = distance[i];
			}
		}
		if (pso.tol > std::abs(f(min_cost) - opt) || rmax < pso.tol)
		{
			std::cout << "Found solution at iteration: " << iter << "." << '\n';
			break;
		}
		auto r = pso.generate_r();
		pso.velocity_update(r, iter);
		for (auto i = 0; i < pso.npop; ++i)
		{
			if (f(pso.particles[i]) < f(pso.personal_best[i]))
			{
				for (auto j = 0; j < pso.ndv; j++)
				{
					pso.personal_best[i][j] = pso.particles[i][j];
				}
			}
		}
		for (auto k = 0; k < pso.nneigh; ++k)
		{
			for (auto l = 0; l < pso.neighbourhoods[k].size(); ++l)
			{
				auto i = pso.neighbourhoods[k][l];
				if (f(pso.personal_best[i]) < f(pso.local_best[k]))
				{
					for (auto j = 0; j < pso.ndv; j++)
					{
						pso.local_best[k][j] = pso.personal_best[i][j];
					}
				}
			}
		}
		min_cost = find_min_local_best(pso.local_best, f);
	}
	return min_cost;
}