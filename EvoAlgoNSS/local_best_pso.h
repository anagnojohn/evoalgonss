#pragma once

#include "dependencies.h"

template<typename T>
class LocalBestPSO
{
public:
	LocalBestPSO(std::vector<std::vector<T>> i_particles, const T& i_tol, const size_t& i_iter_max, const T& i_c1, const T& i_c2, const T& i_alpha, const T& i_w,const std::vector<T>& i_stdev)
		: particles(i_particles), c1(i_c1), c2(i_c2), alpha(i_alpha), tol(i_tol), iter_max(i_iter_max), stdev(i_stdev), npop(particles.size()), ndv(particles[0].size())
	{
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		init_epsilon(particles, stdev);
	}
	T euclid_distance(const std::vector<T>& x, const std::vector<T>& y);
	std::vector<std::vector<T>> generate_r();
	const T c1; //1;
	const T c2; //1;
	const T alpha; //2;
	T w; //0.9;
	std::vector<std::vector<T>> particles;
	const size_t iter_max;
	const T tol;
	const std::vector<T> stdev;
	size_t npop;
	const size_t ndv;
private:
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
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
T LocalBestPSO<T>::euclid_distance (const std::vector<T>& x, const std::vector<T>& y)
{
	T sum = 0;
	for (auto i = 0; i < x.size(); ++i)
	{
		sum = sum + std::pow(x[i] - y[i], 2);
	}
	return std::sqrt(sum);
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
	std::vector<T> vmax(pso.ndv);
	for (auto& p : vmax)
	{
		p = 1000;
	}
	T rmax;
	size_t nneigh = static_cast<size_t>(std::ceil(pso.npop / 5));
	std::vector<T> distance;
	std::vector<std::vector<T>> personal_best;
	std::vector<std::vector<T>> local_best(nneigh);
	std::vector<std::vector<T>> velocity(pso.npop, std::vector<T>(pso.ndv));
	
	std::vector<std::vector<size_t>> neighbourhoods(nneigh);
	std::vector<size_t> indices;
	std::vector<T> min_cost(pso.ndv);
	for (auto i = 0; i < pso.npop; ++i)
	{
		indices.push_back(i);
	}
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	size_t counter = 0;
	for (auto k = 0; k < nneigh; ++k)
	{
		while (counter < 5)
		{
			neighbourhoods[k].push_back(indices[k * 5 + counter]);
			counter = counter + 1;
		}
	}
	// Change when developing final version
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto j = 0; j < pso.ndv; ++j)
		{
			local_best[k].push_back(pso.particles[0][j]);
		}
	}
	for (auto& p : velocity)
	{
		for (auto& n : p)
		{
			n = 0.0;
		}
	}
	for (const auto& p : pso.particles)
	{
		personal_best.push_back(p);
	}
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			if (f(personal_best[i]) < f(local_best[k]))
			{
				for (auto j = 0; j < pso.ndv; j++)
				{
					local_best[k][j] = personal_best[i][j];
				}
			}
		}
	}
	min_cost = find_min_local_best(local_best, f);
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
		for (auto k = 0; k < nneigh; ++k)
		{
			for (auto l = 0; l < neighbourhoods[k].size(); ++l)
			{
				auto i = neighbourhoods[k][l];
				for (auto j = 0; j < pso.ndv; ++j)
				{
					velocity[i][j] = pso.w * velocity[i][j] + pso.c1 * r[0][j] * (personal_best[i][j]
						- pso.particles[i][j]) + pso.c2 * r[1][j] * (local_best[k][j] - pso.particles[i][j]);
					if (velocity[i][j] > vmax[j])
					{
						velocity[i][j] = vmax[j];
					}
					pso.particles[i][j] = pso.particles[i][j] + velocity[i][j];
				}
			}
		}
		for (auto i = 0; i < pso.npop; ++i)
		{
			if (f(pso.particles[i]) < f(personal_best[i]))
			{
				for (auto j = 0; j < pso.ndv; j++)
				{
					personal_best[i][j] = pso.particles[i][j];
				}
			}
		}
		for (auto k = 0; k < nneigh; ++k)
		{
			for (auto l = 0; l < neighbourhoods[k].size(); ++l)
			{
				auto i = neighbourhoods[k][l];
				if (f(personal_best[i]) < f(local_best[k]))
				{
					for (auto j = 0; j < pso.ndv; j++)
					{
						local_best[k][j] = personal_best[i][j];
					}
				}
			}
		}
		min_cost = find_min_local_best(local_best, f);
		for (auto& p : vmax)
		{
			p = (1 - std::pow(iter / pso.iter_max, pso.alpha)) * p;
		}
		pso.w = ((pso.w - 0.4) * (pso.iter_max - iter)) / (pso.iter_max + 0.4);
	}
	return min_cost;
}