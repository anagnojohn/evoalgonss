#pragma once

#include "dependencies.h"

template<typename T, typename F>
class GeneticAlgo
{
public:
	GeneticAlgo(std::vector<std::vector<T>> i_individuals, const T& i_tol, const T& i_opt, const size_t& i_gmax, const T& i_x_rate, const T& i_pi, const std::vector<T>& i_stdev)
		: individuals(i_individuals), x_rate(i_x_rate), pi(i_pi), tol(i_tol), opt(i_opt), gmax(i_gmax), stdev(i_stdev), npop(individuals.size()), ndv(individuals[0].size())
	{
		boost::math::beta_distribution<T> i_dist(1, 6);
		dist = i_dist;
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	std::vector<T> genetic_algo(F f);
private:
	std::vector<std::vector<T>> individuals;
	const T x_rate;// = 0.4;
	const T pi;// = 0.35;
	const T tol;
	const T opt;
	std::vector<T> stdev;
	const size_t gmax;
	size_t npop;
	const size_t ndv;
	std::random_device generator;
	boost::math::beta_distribution<T> dist;
	std::uniform_real_distribution<T> distribution;
	std::vector<T> blend(std::vector<T> r, std::vector<T> s);
	std::vector<std::vector<T>> selection(F f);
	void mutation();
};

template<typename T, typename F>
std::vector<T> GeneticAlgo<T,F>::blend(std::vector<T> r, std::vector<T> s)
{
	std::vector<T> offspring(ndv);
	std::vector<T> psi(ndv);
	for (auto i = 0; i < ndv; ++i)
	{
		psi[i] = distribution(generator);
	}
	for (auto i = 0; i < ndv; ++i)
	{
		offspring[i] = psi[i] * r[i] + (1 - psi[i])*s[i];
	}
	return offspring;
}

template<typename T, typename F>
std::vector<std::vector<T>> GeneticAlgo<T, F>::selection(F f)
{
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

template<typename T, typename F>
void GeneticAlgo<T, F>::mutation()
{
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
std::vector<T> GeneticAlgo<T, F>::genetic_algo(F f)
{
	T best_cost;
	size_t nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	init_epsilon(individuals, stdev);
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
		std::vector<std::vector<T>> offspring = selection(f);
		individuals.erase(individuals.begin() + nkeep, individuals.begin() + npop);
		npop = individuals.size();
		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
		for (auto i = 0; i < offspring.size(); ++i)
		{
			individuals.push_back(offspring[i]);
		}
		mutation();
		for (auto j = 0; j < ndv; ++j)
		{
			stdev[j] = stdev[j] + 0.02 * stdev[j];
		}
	}
	std::sort(individuals.begin(), individuals.end(), comparator);
	return individuals[0];
}