#pragma once

#include "dependencies.h"

template<typename T>
std::tuple<T, T, T> pso()
{
	auto comparator = [&](const std::tuple<T, T, T>& l, const std::tuple<T, T, T>& r)
	{
		return std::get<2>(l) < std::get<2>(r);
	};
	std::default_random_engine generator;
	std::uniform_real_distribution<T> distribution(0.0, 1.0);
	T inf = std::numeric_limits<T>::infinity();
	T c1 = 1;
	T c2 = 1;
	T npop = 1000;
	T iter_max = 1000;
	T tol = 0.001;
	std::vector<std::tuple<T, T, T>> local_best;
	std::vector<std::tuple<T, T>> velocity;
	std::tuple<T, T, T> global_best = { inf, inf, inf };
	T dc1_max = 10;
	T dc1_min = 0;
	T dc2_max = 10;
	T dc2_min = 0;
	std::vector<T> dist;
	T rmax = 0.0;
	T opt = 1.0;
	auto ndv = 2;
	auto euclid_distance = [&](std::tuple<T, T, T> x, std::tuple<T, T, T> y)
	{ return std::sqrt(std::pow(std::get<0>(x) - std::get<0>(y), 2) + std::pow(std::get<1>(x) - std::get<1>(y), 2)); };
	auto f = [&](std::tuple<T, T, T> dv) { return std::pow(std::get<0>(dv), 2) + std::pow(std::get<1>(dv), 2) + 1; };
	std::vector<std::tuple<T, T, T>> particles;
	for (auto i = 0; i < npop; ++i)
	{
		T x = dc1_min + (distribution(generator)) * (dc1_max - dc1_min);
		T y = dc2_min + (distribution(generator)) * (dc2_max - dc2_min);
		particles.push_back({ x,y,0.0 });
		velocity.push_back({ 0.0, 0.0 });
	}
	for (auto& p : particles)
	{
		std::get<2>(p) = f(p);
	}
	for (const auto& p : particles)
	{
		local_best.push_back(p);
	}
	for (auto i = 0; i < npop; ++i)
	{
		std::get<2>(local_best[i]) = f(local_best[i]);
		if (std::get<2>(local_best[i]) < std::get<2>(global_best))
		{
			std::get<0>(global_best) = std::get<0>(local_best[i]);
			std::get<1>(global_best) = std::get<1>(local_best[i]);
			std::get<2>(global_best) = std::get<2>(local_best[i]);
		}
	}
	std::cout << "particles:" << std::endl;
	for (const auto& p : particles)
	{
		std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
	}
	for (auto iter = 0; iter < iter_max; ++iter)
	{
		if (tol > std::abs(std::get<2>(global_best) - opt))
		{
			std::cout << "Found solution at iteration: " << iter << "." << std::endl;
			break;
		}
		for (auto i = 0; i < npop; ++i)
		{
			dist.push_back(euclid_distance(particles[i], global_best));
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
		T r11 = distribution(generator);
		T r12 = distribution(generator);
		T r21 = distribution(generator);
		T r22 = distribution(generator);
		for (auto i = 0; i < npop; ++i)
		{	
			std::get<0>(velocity[i]) = std::get<0>(velocity[i]) + c1 * r11 * (std::get<0>(local_best[i])
				- std::get<0>(particles[i])) + c2 * r21 * (std::get<0>(global_best) - std::get<0>(particles[i]));
			std::get<1>(velocity[i]) = std::get<1>(velocity[i]) + c1 * r12 * (std::get<1>(local_best[i])
				- std::get<1>(particles[i])) + c2 * r22 * (std::get<1>(global_best) - std::get<1>(particles[i]));
			std::get<0>(particles[i]) = std::get<0>(particles[i]) + std::get<0>(velocity[i]);
			std::get<1>(particles[i]) = std::get<1>(particles[i]) + std::get<1>(velocity[i]);
			std::get<2>(particles[i]) = f(particles[i]);
		}
		for (auto i = 0; i < npop; ++i)
		{
			if (std::get<2>(particles[i]) < std::get<2>(local_best[i]))
			{
				std::get<0>(local_best[i]) = std::get<0>(particles[i]);
				std::get<1>(local_best[i]) = std::get<1>(particles[i]);
				std::get<2>(local_best[i]) = std::get<2>(particles[i]);
			}
			if (std::get<2>(local_best[i]) < std::get<2>(global_best))
			{
				std::get<0>(global_best) = std::get<0>(local_best[i]);
				std::get<1>(global_best) = std::get<1>(local_best[i]);
				std::get<2>(global_best) = std::get<2>(local_best[i]);
			}
		}
	}
	std::sort(particles.begin(), particles.end(), comparator);
	std::cout << "Sorted Costs:" << std::endl;
	for (const auto& p : particles)
	{
		std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
	}
	return global_best;
}
