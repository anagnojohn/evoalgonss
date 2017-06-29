#pragma once

#include "dependencies.h"

template<typename T, typename F, typename U>
std::vector<std::vector<T>> selection(const std::vector<std::vector<T>>& individuals, F f, size_t nkeep, U comparator)
{
	auto npop = individuals.size();
	auto ndv = individuals[0].size();
	std::default_random_engine generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	boost::math::beta_distribution<> dist(1, 6);
	T optimal_cost = f(individuals[0]);
	T sum_cost = 0.0;
	for (auto i = 0; i < nkeep; ++i)
	{
		sum_cost = sum_cost + f(individuals[i]) - f(individuals[nkeep]);
	}
	std::vector<T> prob;
	for (auto i = 0; i < nkeep; ++i)
	{
		prob.push_back(std::abs((f(individuals[i]) - f(individuals[nkeep])) / sum_cost));
	}
	std::vector<T> cum_prob;
	cum_prob.push_back(prob[0]);
	for (auto i = 1; i < nkeep; ++i)
	{
		cum_prob.push_back(cum_prob[i - 1] + prob[i]);
	}
	std::vector<std::vector<T>> parents;
	for (auto i = 0; i < nkeep; ++i)
	{
		T xi = quantile(dist, distribution(generator));
		T r = xi;
		if (r < cum_prob[i])
		{
			parents.push_back(individuals[i]);
		}
	}
	std::sort(parents.begin(), parents.end(), comparator);
	if (parents.size() % 2 != 0)
	{
		parents.erase(parents.begin() + parents.size() - 1, parents.end());
	}
	return parents;
}

template<typename T>
std::vector<T> blend(std::vector<T> psi, std::vector<T> r, std::vector<T> s)
{
	std::vector<T> offspring(psi.size());
	for (auto i = 0; i < psi.size(); ++i)
	{
		offspring[i] = psi[i] * r[i] + (1 - psi[i])*s[i];
	}
	return offspring;
}

template<typename T>
void mutation(std::vector<std::vector<T>>& individuals, const std::vector<T>& epsilon)
{
	std::default_random_engine generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	auto npop = individuals.size();
	auto ndv = individuals[0].size();
	T pi = 0.35;
	for (auto i = 1; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			T r = distribution(generator);
			if (pi < r)
			{
				individuals[i][j] = individuals[i][j] + epsilon[j];
			}

		}
	}
}

template<typename T, typename F>
std::vector<T> genetic_algo(std::vector<std::vector<T>> individuals, F f)
{
	size_t npop = individuals.size();
	size_t ndv = individuals[0].size();
	T gmax = 5000;
	T x_rate = 0.5;
	size_t nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	T tol = 0.00001;
	T opt = 1.0;

	std::default_random_engine generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	boost::math::beta_distribution<> dist(1, 6);
	//auto epsilon = init_epsilon(individuals);
	std::vector<T> epsilon;
	for (auto j = 0; j < ndv; ++j)
	{
		epsilon.push_back(0.5);
	}
	//std::cout << "Individuals:" << std::endl;
	//for (const auto& p : individuals)
	//{
	//	std::cout << p << " " << f(p) << std::endl;
	//}
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			individuals[i][j] = individuals[i][j] + epsilon[j];
		}
	}
	for (auto g = 0; g < gmax; ++g)
	{
		if (tol > std::abs(f(individuals[0]) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << std::endl;
			break;
		}
		std::sort(individuals.begin(), individuals.end(), comparator);

		//std::cout << parents.size();
		std::vector<T> psi(ndv);
		for (auto i = 0; i < ndv; ++i)
		{
			psi[i] = distribution(generator);
		}
		std::vector<std::vector<T>> parents = selection(individuals, f, nkeep, comparator);
		individuals.erase(individuals.begin() + nkeep, individuals.begin() + npop);
		for (auto i = 0; i < parents.size(); i += 2)
		{
			individuals.push_back(blend(psi, parents[i], parents[i + 1]));
		}
		npop = individuals.size();
		if (individuals.size() < 3)
		{
			break;
		}

		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
		mutation(individuals, epsilon);
		for (auto j = 0; j < ndv; ++j)
		{
			epsilon[j] = epsilon[j] + 0.02 * epsilon[j];
		}
	}
	std::sort(individuals.begin(), individuals.end(), comparator);
	std::cout << "Sorted Costs:" << std::endl;
	for (const auto& p : individuals)
	{
		std::cout << p << " " << f(p) << std::endl;
	}
	return individuals[0];
}