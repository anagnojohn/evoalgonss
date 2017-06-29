#pragma once

#include "dependencies.h"

template<typename T>
T euclid_distance(const std::vector<T>& x, const std::vector<T>& y)
{
	T sum = 0;
	for (auto i = 0; i < x.size(); ++i)
	{
		sum = sum + std::pow(x[i] - y[i], 2);
	}
	return std::sqrt(sum);
};

template<typename T, typename F>
std::vector<T> gbest_pso(std::vector<std::vector<T>> particles, F f)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<T> distribution(0.0, 1.0);
	T inf = std::numeric_limits<T>::infinity();
	T c1 = 1;
	T c2 = 1;
	size_t npop = particles.size();
	size_t ndv = particles[0].size();
	T iter_max = 1000;
	T tol = 0.001;
	T rmax = 0.0;
	T opt = 1.0;
	T vmax = 0.5;
	std::vector<T> dist;
	std::vector<std::vector<T>> local_best;
	std::vector<std::vector<T>> velocity(npop, std::vector<T>(ndv));
	std::vector<T> global_best = { inf, inf };
	std::vector<std::vector<T>> r(2, std::vector<T>(ndv));
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	// Change when developing final version
	for (auto& p : velocity)
	{
		for (auto& n : p)
		{
			n = 0.0;
		}
	}
	for (const auto& p : particles)
	{
		local_best.push_back(p);
	}
	for (auto i = 0; i < npop; ++i)
	{
		if (f(local_best[i]) < f(global_best))
		{
			for (auto j = 0; j < ndv; j++)
			{
				global_best[j] = local_best[i][j];
			}
		}
	}
	for (const auto& p : particles)
	{
		//	std::cout << p << " " << f(p) << std::endl;
	}
	for (auto iter = 0; iter < iter_max; ++iter)
	{
		if (tol > std::abs(f(global_best) - opt))
		{
			std::cout << "Found solution at iteration: " << iter << "." << std::endl;
			break;
		}
		for (auto i = 0; i < npop; ++i)
		{
			dist.push_back(euclid_distance<T>(particles[i], global_best));
		}
		rmax = dist[0];
		for (auto i = 0; i < npop; ++i)
		{
			if (rmax < dist[i])
			{
				rmax = dist[i];
			}
		}
		if (rmax < tol)
		{
			std::cout << "Found solution at iteration: " << iter << "." << std::endl;
			break;
		}
		for (auto i = 0; i < 2; ++i)
		{
			for (auto j = 0; j < ndv; ++j)
			{
				r[i][j] = (distribution(generator));
			}
		}
		for (auto i = 0; i < npop; ++i)
		{
			for (auto j = 0; j < ndv; ++j)
			{
				velocity[i][j] = velocity[i][j] + c1 * r[0][j] * (local_best[i][j]
					- particles[i][j]) + c2 * r[1][j] * (global_best[j] - particles[i][j]);
				particles[i][j] = particles[i][j] + velocity[i][j];
			}
		}
		for (auto i = 0; i < npop; ++i)
		{
			if (f(particles[i]) < f(local_best[i]))
			{
				for (auto j = 0; j < ndv; j++)
				{
					local_best[i][j] = particles[i][j];
				}
			}
			if (f(local_best[i]) < f(global_best))
			{
				for (auto j = 0; j < ndv; j++)
				{
					global_best[j] = local_best[i][j];
				}
			}
		}
	}
	std::sort(particles.begin(), particles.end(), comparator);
	std::cout << "Sorted Costs:" << std::endl;
	for (const auto& p : particles)
	{
		//	std::cout << p << " " << f(p) << std::endl;
	}
	return global_best;
}