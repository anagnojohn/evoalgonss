#pragma once

#include <vector>
#include <random>
#include <boost/math/distributions.hpp>
#include <tuple>
#include <chrono>
#include <ctime>
#include "ealgorithm_base.h"

// Genetic Algorithms Class
template<typename T>
class GeneticAlgo : public EA_base<T>
{
public:
	GeneticAlgo(const T& i_x_rate, const T& i_pi, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: x_rate{ i_x_rate }, pi{ i_pi }, EA_base{ i_npop, i_tol, i_iter_max }
	{
		assert(x_rate > 0 && x_rate <= 1);
		assert(pi > 0 && pi <= 1);
		assert(tol > 0);
		assert(iter_max > 0);
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
std::tuple<std::vector<T>, T, size_t, double> solve(F f, const T& opt, const GeneticAlgo<T>& ga, const EAparams<T>& ea)
{
	std::cout << "Genetic Algorithm used as solver" << "\n";
	// Find the minimum cost individual of the fitness function for the population
	const auto& tol = ga.get_tol();
	const auto& iter_max = ga.get_iter_max();
	auto npop = ea.get_npop();
	const auto& ndv = ea.get_ndv();
	auto stdev = ea.get_stdev();
	auto individuals = ea.get_individuals();
	std::vector<T> min_cost = individuals[0];
	const auto& x_rate = ga.get_x_rate();
	size_t nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	T fitness_cost = f(individuals[0]);
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	size_t last_iter = 0;
	//  Time the computation
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	std::sort(individuals.begin(), individuals.end(), comparator);
	for (auto iter = 0; iter < iter_max; ++iter)
	{
		if (tol > std::abs(fitness_cost - opt))
		{
			break;
		}
		if (individuals.size() < 3)
		{
			last_iter = iter;
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
		std::sort(individuals.begin(), individuals.end(), comparator);
		fitness_cost = f(individuals[0]);
		last_iter = iter;
	}
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	T timer = elapsed_seconds.count();
	return { individuals[0], fitness_cost, last_iter, timer };
}