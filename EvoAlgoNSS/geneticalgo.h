#pragma once

#include <boost/math/distributions.hpp>
#include "ealgorithm_base.h"

template<typename T>
struct GAstruct : EAstruct<T>
{
public:
	GAstruct(const T& i_x_rate, const T& i_pi, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: x_rate{ i_x_rate }, pi{ i_pi }, EAstruct{ i_npop, i_tol, i_iter_max }
	{
		assert(x_rate > 0 && x_rate <= 1);
		assert(pi > 0 && pi <= 1);
	}
	// Natural Selection rate
	const T x_rate;
	// Probability of mutating
	const T pi;
};

// Genetic Algorithms Class
template<typename T, typename F>
class Solver<T, F, GAstruct<T>> : public Solver<T, F, EAstruct<T>>
{
public:
	Solver(const GAstruct<T>& ga, const Population<T>& popul) : Solver < T, F, EAstruct<T>> { { ga.npop, ga.tol, ga.iter_max }, popul}, x_rate{ ga.x_rate }, pi{ ga.pi }, stdev{ popul.stdev }
	{
		boost::math::beta_distribution<T> i_dist(1, 6);
		dist = i_dist;
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	}
private:
	// The standard deviation of variables is not constant in GA
	std::vector<T> stdev;
	std::vector<size_t> indices;
	// Natural Selection rate
	T x_rate;
	// Probability of mutating
	T pi;
	// Number of individuals to be kept
	size_t nkeep;
	std::random_device generator;
	boost::math::beta_distribution<T> dist;
	std::uniform_real_distribution<T> distribution;
	std::string get_type_of_solver() { return "Genetic Algorithms"; };
	// Crossover method
	std::vector<T> blend(std::vector<T> r, std::vector<T> s);
	// Selection method
	std::vector<std::vector<T>> selection();
	void mutation();
	void run_algo(F f, const T& opt);
};

template<typename T, typename F>
std::vector<T> Solver<T, F, GAstruct<T>>::blend(std::vector<T> r, std::vector<T> s)
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
std::vector<std::vector<T>> Solver<T, F, GAstruct<T>>::selection()
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

template<typename T, typename F>
void Solver<T, F, GAstruct<T>>::mutation()
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
void Solver<T, F, GAstruct<T>>::run_algo(F f, const T& opt)
{
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	std::sort(individuals.begin(), individuals.end(), comparator);
	for (auto iter = 0; iter < iter_max; ++iter)
	{
		if (tol > std::abs(fitness_cost - opt))
		{
			last_iter = iter;
			break;
		}
		if (individuals.size() < 3)
		{
			last_iter = iter;
			break;
		}
		std::vector<std::vector<T>> offspring = selection();
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
		std::sort(individuals.begin(), individuals.end(), comparator);
		fitness_cost = f(individuals[0]);
		last_iter = iter;
	}
}