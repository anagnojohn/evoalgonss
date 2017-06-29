#pragma once

#include "dependencies.h"

template<typename T>
std::tuple<T, T, T> genetic_algo()
{
	T beta = 0.25;
	T npop = 6;
	T crate = 0.8;
	T mrate = 0.2;
	T gmax = 1000;
	T nkeep = npop - 2;
	T tol = 0.00001;
	//auto print = [](const int& n) { std::cout << " " << n; };
	//std::vector<T> cost = (individuals.begin(), individuals.end(), [](T a) { return std::toupper(c); });
	//std::cout << "before:";
	// or std::cout << std::get<0>(p) << ", " << std::get<1>(p) << std::endl;
	//std::for_each(nums.begin(), nums.end(), print);
	//std::cout << '\n';
	auto f = [&](std::tuple<T, T, T> dv) { return std::pow(std::get<0>(dv), 2) + std::pow(std::get<1>(dv), 2) + 1; };
	std::vector<std::tuple<T, T, T>> individuals;
	for (auto i = 0; i < npop; ++i)
	{
		T x = rand() % 10;
		T y = rand() % 10;
		individuals.push_back({ x,y,0.0 });
	}
	for (auto& p : individuals)
	{
		std::get<2>(p) = f(p);
	}
	std::cout << "Individuals:" << std::endl;
	for (const auto& p : individuals)
	{
		std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
	}
	auto comparator = [&](const std::tuple<T, T, T>& l, const std::tuple<T, T, T>& r)
	{
		return std::get<2>(l) < std::get<2>(r);
	};
	for (auto g = 0; g < gmax; ++g)
	{
		if (tol > std::get<2>(individuals[0]))
		{
			std::cout << "Found solution at iteration: " << g << "." << std::endl;
			break;
		}
		//std::vector<T> costs = std::for_each(individuals.begin(), individuals.end(), f);
		std::sort(individuals.begin(), individuals.end(), comparator);
		//std::cout << "Sorted Costs:" << std::endl;
		//for (const auto& p : individuals)
		//{
		//	std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
		//}
		for (auto i = 0; i < nkeep; ++i)
		{
			std::get<2>(individuals[i]) = std::get<2>(individuals[i]) - std::get<2>(individuals[nkeep + 1]);
		}
		//std::cout << "Costs Weighting:" << std::endl;
		//for (const auto& p : individuals)
		//{
		//	std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
		//}
		T sum = 0.0;
		for (auto i = 0; i < nkeep; ++i)
		{
			sum = sum + std::get<2>(individuals[i]);
		}
		std::vector<T> prob;
		for (auto i = 0; i < nkeep; ++i)
		{
			prob.push_back(std::get<2>(individuals[i]) / sum);
		}
		//sum = 0.0;
		//std::cout << "Probabilities" << std::endl;
		//for (const auto& p : prob)
		//{
		//	std::cout << p << std::endl;
		//	sum = sum + p;
		//}
		std::vector<T> cum_prob;
		cum_prob.push_back(prob[0]);
		for (auto i = 1; i < nkeep; ++i)
		{
			cum_prob.push_back(cum_prob[i - 1] + prob[i]);
		}
		//std::cout << "Cumulative Probabilities" << std::endl;
		//for (const auto& p : cum_prob)
		//{
		//	std::cout << p << std::endl;
		//}
		std::vector<std::tuple<T, T, T>> parents;
		auto np = 0;
		for (auto i = 0; i < nkeep; ++i)
		{
			T r = rand() % 1;
			if (r < cum_prob[i])
			{
				if (np > 1)
				{
					break;
				}
				parents.push_back(individuals[i]);
				np = np + 1;
			}
		}
		//std::cout << "Parents:" << std::endl;
		//for (const auto& p : parents)
		//{
		//	std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
		//}
		for (auto i = nkeep; i < npop; i += 2)
		{
			std::get<0>(individuals[i]) = beta*std::get<0>(parents[0]) + (1 - beta)*std::get<0>(parents[1]);
			std::get<1>(individuals[i]) = beta*std::get<1>(parents[0]) + (1 - beta)*std::get<1>(parents[1]);
			std::get<0>(individuals[i + 1]) = (1 - beta)*std::get<0>(parents[0]) + beta*std::get<0>(parents[1]);
			std::get<1>(individuals[i + 1]) = (1 - beta)*std::get<1>(parents[0]) + beta*std::get<1>(parents[1]);
		}
		for (auto& p : individuals)
		{
			std::get<2>(p) = f(p);
		}
		//std::cout << "Population after crossover: " << std::endl;
		//for (const auto& p : individuals)
		//{
		//	std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
		//}
		for (auto i = 1; i < npop; i += 3)
		{
			T r = rand() % 10;
			std::get<0>(individuals[i]) = r;
			r = rand() % 10;
			std::get<1>(individuals[i]) = r;
			std::get<2>(individuals[i]) = f(individuals[i]);
		}
		//std::cout << "Population after mutation: " << std::endl;
		//for (const auto& p : individuals)
		//{
		//	std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
		//}
	}
	std::sort(individuals.begin(), individuals.end(), comparator);
	std::cout << "Sorted Costs:" << std::endl;
	for (const auto& p : individuals)
	{
		std::cout << std::get<0>(p) << ", " << std::get<1>(p) << ", " << std::get<2>(p) << std::endl;
	}
	return individuals[0];
}