#pragma once

#include "dependencies.h"

template<typename T, typename F>
std::vector<T> find_min_local_best(const std::vector<std::vector<T>>& local_best, F f)
{
	T inf = std::numeric_limits<T>::infinity();
	auto nneigh = local_best.size();
	auto ndv = local_best[0].size();
	std::vector<T> min_cost(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		min_cost[j] = inf;
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

template<typename T>
T euclid_distance (const std::vector<T>& x, const std::vector<T>& y)
{
	T sum = 0;
	for (auto i = 0; i < x.size(); ++i)
	{
		sum = sum + std::pow(x[i] - y[i], 2);
	}
	return std::sqrt(sum);
};

template<typename T, typename F>
std::vector<T> lbest_pso(std::vector<std::vector<T>> particles, F f, T tol, T opt, size_t iter_max)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<T> distribution(0.0, 1.0);
	T inf = std::numeric_limits<T>::infinity();
	T c1 = 1;
	T c2 = 1;
	size_t npop = particles.size();
	size_t ndv = particles[0].size();
	T rmax = 0.0;
	T alpha = 2;
	T w = 0.9;
	std::vector<T> vmax(ndv);
	for (auto& p : vmax)
	{
		p = 1000;
	}
	size_t nneigh = static_cast<size_t>(std::ceil(npop / 5));
	std::vector<T> dist;
	std::vector<std::vector<T>> personal_best;
	std::vector<std::vector<T>> local_best(nneigh);
	std::vector<std::vector<T>> velocity(npop, std::vector<T>(ndv));
	std::vector<std::vector<T>> r(2, std::vector<T>(ndv));
	std::vector<std::vector<size_t>> neighbourhoods(nneigh);
	std::vector<size_t> indices;
	std::vector<T> min_cost(ndv);
	for (auto i = 0; i < npop; ++i)
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
		for (auto j = 0; j < ndv; ++j)
		{
			local_best[k].push_back(inf);
		}
	}
	for (auto& p : velocity)
	{
		for (auto& n : p)
		{
			n = 0.0;
		}
	}
	for (const auto& p : particles)
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
				for (auto j = 0; j < ndv; j++)
				{
					local_best[k][j] = personal_best[i][j];
				}
			}
		}
	}
	for (const auto& p : particles)
	{
	//	std::cout << p << " " << f(p) << '\n';
	}
	min_cost = find_min_local_best(local_best, f);
	for (auto iter = 0; iter < iter_max; ++iter)
	{	
		for (auto i = 0; i < npop; ++i)
		{
			dist.push_back(euclid_distance<T>(particles[i], min_cost));
		}
		rmax = dist[0];
		for (auto i = 0; i < npop; ++i)
		{
			if (rmax < dist[i])
			{
				rmax = dist[i];
			}
		}
		if (tol > std::abs(f(min_cost) - opt) || rmax < tol)
		{
			std::cout << "Found solution at iteration: " << iter << "." << '\n';
			break;
		}
		for (auto i = 0; i < 2; ++i)
		{
			for (auto j = 0; j < ndv; ++j)
			{
				r[i][j] = (distribution(generator));
			}
		}
		for (auto k = 0; k < nneigh; ++k)
		{
			for (auto l = 0; l < neighbourhoods[k].size(); ++l)
			{
				auto i = neighbourhoods[k][l];
				for (auto j = 0; j < ndv; ++j)
				{
					velocity[i][j] = w * velocity[i][j] + c1 * r[0][j] * (personal_best[i][j]
						- particles[i][j]) + c2 * r[1][j] * (local_best[k][j] - particles[i][j]);
					if (velocity[i][j] > vmax[j])
					{
						velocity[i][j] = vmax[j];
					}
					particles[i][j] = particles[i][j] + velocity[i][j];
				}
			}
		}
		for (auto i = 0; i < npop; ++i)
		{
			if (f(particles[i]) < f(personal_best[i]))
			{
				for (auto j = 0; j < ndv; j++)
				{
					personal_best[i][j] = particles[i][j];
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
					for (auto j = 0; j < ndv; j++)
					{
						local_best[k][j] = personal_best[i][j];
					}
				}
			}
		}
		min_cost = find_min_local_best(local_best, f);
		std::cout << f(min_cost) << '\n';
		for (auto& p : vmax)
		{
			p = (1 - std::pow(iter / iter_max, alpha)) * p;
		}
		w = ((w - 0.4) * (iter_max - iter)) / (iter_max + 0.4);
	}
	//std::sort(particles.begin(), particles.end(), comparator);
	//std::cout << "Sorted Costs:" << '\n';
	for (const auto& p : particles)
	{
	//	std::cout << p << " " << f(p) << '\n';
	}
	//std::sort(local_best.begin(), local_best.end(), comparator);
	for (const auto& p : local_best)
	{
	//	std::cout << p << " " << f(p) << '\n';
	}
	return min_cost;
}