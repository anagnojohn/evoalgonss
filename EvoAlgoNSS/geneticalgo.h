#pragma once

#include "dependencies.h"

template<typename T>
std::vector<T> blend(std::vector<T> r, std::vector<T> s)
{
	auto ndv = r.size();
	std::default_random_engine generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	std::vector<T> psi(ndv);
	for (auto i = 0; i < ndv; ++i)
	{
		psi[i] = distribution(generator);
	}
	std::vector<T> offspring(psi.size());
	for (auto i = 0; i < ndv; ++i)
	{
		offspring[i] = psi[i]*r[i] + (1 - psi[i])*s[i];
	}
	return offspring;
}

template<typename T, typename F>
std::vector<std::vector<T>> selection(const std::vector<std::vector<T>>& individuals, F f)
{
	auto npop = individuals.size();
	auto ndv = individuals[0].size();
	std::random_device generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	boost::math::beta_distribution<> dist(1, 6);
	T optimal_cost = f(individuals[0]);
	std::vector<std::vector<T>> offspring;
	for (auto i = 0; i < npop; ++i)
	{
		T xi = quantile(dist, distribution(generator));
		size_t r = static_cast<size_t>(std::round(npop * xi));
		xi = quantile(dist, distribution(generator));
		size_t s = static_cast<size_t>(std::round(npop * xi));
		{
			offspring.push_back(blend(individuals[r], individuals[s]));
		}
	}
	return offspring;
}

template<typename T>
void mutation(std::vector<std::vector<T>>& individuals, const std::vector<T>& stdev)
{
	std::random_device generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	auto npop = individuals.size();
	auto ndv = individuals[0].size();
	T pi = 0.35;
	T epsilon;
	for (auto i = 1; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			T r = distribution(generator);
			if (pi < r)
			{
				std::normal_distribution<> ndistribution(0, stdev[j]);
				epsilon = ndistribution(generator);
				individuals[i][j] = individuals[i][j] + epsilon;
			}

		}
	}
}

template<typename T, typename F>
std::vector<T> genetic_algo(std::vector<std::vector<T>> individuals, F f, const T& tol, const T& opt, const size_t& gmax, std::vector<T>& stdev)
{
	size_t npop = individuals.size();
	const size_t& ndv = individuals[0].size();
	T x_rate = 0.4;
	size_t nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	T best_cost;
	std::random_device generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	boost::math::beta_distribution<> dist(1, 6);
	init_epsilon(individuals, stdev);
	//std::cout << "Individuals:" << '\n';
	//for (const auto& p : individuals)
	//{
	//	std::cout << p << " " << f(p) << '\n';
	//}
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	init_epsilon(individuals, stdev);
	for (auto g = 0; g < gmax; ++g)
	{
		std::sort(individuals.begin(), individuals.end(), comparator);
		best_cost = f(individuals[0]);
		if (tol > std::abs(f(individuals[0]) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		if (individuals.size() < 3)
		{
			break;
		}
		std::vector<std::vector<T>> offspring = selection(individuals, f);
		individuals.erase(individuals.begin() + nkeep, individuals.begin() + npop);
		npop = individuals.size();
		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
		for (auto i = 0; i < offspring.size(); ++i)
		{
			individuals.push_back(offspring[i]);
		}
		mutation(individuals, stdev);
		for (auto j = 0; j < ndv; ++j)
		{
			stdev[j] = stdev[j] + 0.02 * stdev[j];
		}
	}
	std::sort(individuals.begin(), individuals.end(), comparator);
	//std::cout << "Sorted Costs:" << '\n';
	for (const auto& p : individuals)
	{
		//std::cout << p << " " << f(p) << '\n';
	}
	return individuals[0];
}