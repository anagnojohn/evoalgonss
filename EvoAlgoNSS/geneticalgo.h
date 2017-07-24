#pragma once

#include "ealgorithm_base.h"

// Genetic Algorithms Class
template<typename T>
class GeneticAlgo : EA_base<T>
{
public:
	GeneticAlgo(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max, const T& i_x_rate, const T& i_pi, const std::vector<T>& stdev)
		: x_rate(i_x_rate), pi(i_pi)
	{
		set_solver(decision_variables, npop, tol, iter_max);
		boost::math::beta_distribution<T> i_dist(1, 6);
		dist = i_dist;
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	std::vector<T> blend(std::vector<T> r, std::vector<T> s);
	std::vector<std::vector<T>> selection();
	void mutation();
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
std::vector<std::vector<T>> GeneticAlgo<T>::selection()
{
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
void GeneticAlgo<T>::mutation()
{
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
std::vector<T> solve(F f, const T& opt, GeneticAlgo<T>& ga)
{
	T best_cost;
	size_t nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	for (auto g = 0; g < ga.iter_max; ++g)
	{
		std::sort(ga.individuals.begin(), ga.individuals.end(), comparator);
		best_cost = f(ga.individuals[0]);
		if (tol > std::abs(f(ga.individuals[0]) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		if (individuals.size() < 3)
		{
			break;
		}
		std::vector<std::vector<T>> offspring = ga.selection(f);
		ga.individuals.erase(ga.individuals.begin() + nkeep, ga.individuals.begin() + ga.npop);
		npop = ga.individuals.size();
		nkeep = static_cast<size_t>(std::ceil(ga.npop * ga.x_rate));
		for (auto i = 0; i < offspring.size(); ++i)
		{
			ga.individuals.push_back(offspring[i]);
		}
		ga.mutation();
		for (auto j = 0; j < ga.ndv; ++j)
		{
			ga.stdev[j] = ga.stdev[j] + 0.02 * ga.stdev[j];
		}
	}
	std::sort(ga.individuals.begin(), ga.individuals.end(), comparator);
	return ga.individuals[0];
}