#pragma once

#include <vector>
#include <random>
#include <boost/math/distributions.hpp>
#include "ealgorithm_base.h"

// Genetic Algorithms Class
template<typename T>
class GeneticAlgo
{
public:
	GeneticAlgo(const T& i_x_rate, const T& i_pi)
		: x_rate(i_x_rate), pi(i_pi)
	{
		assert((x_rate > 0 && x_rate <= 1) && (pi > 0 && pi <= 1));
		boost::math::beta_distribution<T> i_dist(1, 6);
		dist = i_dist;
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	// Crossover method
	std::vector<T> blend(std::vector<T> r, std::vector<T> s);
	// Selection method
	std::vector<std::vector<T>> selection(const std::vector<std::vector<T>>& individuals);
	void mutation(std::vector<std::vector<T>>& individuals, const std::vector<T>& stdev);
	// Getters
	T get_x_rate() const { return x_rate; }
	// Setters
	void set_x_rate(const T& x_rate) { assert(x_rate > 0 && x_rate <= 1); this->x_rate = x_rate; };
	void set_pi(const T& pi) { assert(pi > 0 && pi <= 1); this->pi = pi; };
private:
	const T x_rate;// = 0.4;
	const T pi;// = 0.35;
	std::random_device generator;
	boost::math::beta_distribution<T> dist;
	std::uniform_real_distribution<T> distribution;
};

template<typename T>
std::vector<T> GeneticAlgo<T>::blend(std::vector<T> r, std::vector<T> s)
{
	const auto& ndv = r.size();
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

template<typename T>
std::vector<std::vector<T>> GeneticAlgo<T>::selection(const std::vector<std::vector<T>>& individuals)
{
	const auto& npop = individuals.size();
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
void GeneticAlgo<T>::mutation(std::vector<std::vector<T>>& individuals, const std::vector<T>& stdev)
{
	const auto& npop = individuals.size();
	const auto& ndv = individuals[0].size();
	T epsilon;
	for (auto i = 1; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			T r = distribution(generator);
			if (pi < r)
			{
				std::normal_distribution<T> ndistribution(0, stdev[j]);
				epsilon = ndistribution(generator);
				individuals[i][j] = individuals[i][j] + epsilon;
			}

		}
	}
}

template<typename T, typename F>
std::vector<T> solve(F f, const T& opt, GeneticAlgo<T>& ga, const EAparams<T>& ea)
{
	// Find the minimum cost individual of the fitness function for the population
	const auto& tol = ea.get_tol();
	const auto& iter_max = ea.get_iter_max();
	auto npop = ea.get_npop();
	const auto& ndv = ea.get_ndv();
	auto stdev = ea.get_stdev();
	auto individuals = ea.get_individuals();
	std::vector<T> min_cost = individuals[0];
	const auto& x_rate = ga.get_x_rate();
	size_t nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	for (auto g = 0; g < iter_max; ++g)
	{
		std::sort(individuals.begin(), individuals.end(), comparator);
		if (tol > std::abs(f(individuals[0]) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		if (individuals.size() < 3)
		{
			break;
		}
		std::vector<std::vector<T>> offspring = ga.selection(individuals);
		individuals.erase(individuals.begin() + nkeep, individuals.begin() + npop);
		npop = individuals.size();
		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
		for (auto i = 0; i < offspring.size(); ++i)
		{
			individuals.push_back(offspring[i]);
		}
		ga.mutation(individuals, stdev);
		for (auto j = 0; j < ndv; ++j)
		{
			stdev[j] = stdev[j] + 0.02 * stdev[j];
		}
	}
	std::sort(individuals.begin(), individuals.end(), comparator);
	return individuals[0];
}