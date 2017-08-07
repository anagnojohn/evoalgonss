#pragma once

#include <boost/math/distributions.hpp>
#include "ealgorithm_base.h"

template<typename T>
struct GAstruct : EAstruct<T>
{
public:
	GAstruct(const T& i_x_rate, const T& i_pi, const T& i_alpha, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: x_rate{ i_x_rate }, pi{ i_pi }, alpha{ i_alpha }, EAstruct{ i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }
	{
		assert(x_rate > 0 && x_rate <= 1);
		assert(pi > 0 && pi <= 1);
	}
	// Natural Selection rate
	const T x_rate;
	// Probability of mutating
	const T pi;
	// Parameter alpha for Beta distribution
	const T alpha;
};

// Genetic Algorithms (GA) Class
template<typename T>
class Solver<T, GAstruct<T>> : public Solver<T, EAstruct<T>>
{
public:
	Solver(const GAstruct<T>& ga) : Solver < T, EAstruct<T>> { { ga.decision_variables, ga.stdev, ga.npop, ga.tol, ga.iter_max } }, x_rate{ ga.x_rate }, pi{ ga.pi }, alpha{ ga.alpha }
	{
		boost::math::beta_distribution<T> i_dist(1, alpha);
		dist = i_dist;
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	}
	const std::string type = "Genetic Algorithms";
	template<typename F, typename C> void run_algo(F f, C c);
private:
	// Natural Selection rate
	T x_rate;
	// Probability of mutating
	T pi;
	// Parameter alpha for Beta distribution
	T alpha;
	// Number of individuals to be kept
	size_t nkeep;
	std::random_device generator;
	boost::math::beta_distribution<T> dist;
	std::uniform_real_distribution<T> distribution;
	// Crossover method
	std::vector<T> crossover(std::vector<T> r, std::vector<T> s);
	// Selection method
	std::vector<T> selection();
	std::vector<T> mutation(const std::vector<T>& individual);
};

// Crossover step of GA
template<typename T>
std::vector<T> Solver<T, GAstruct<T>>::crossover(std::vector<T> r, std::vector<T> s)
{
	std::vector<T> offspring(ndv);
	std::vector<T> psi(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		psi[j] = distribution(generator);
	}
	for (auto j = 0; j < ndv; ++j)
	{
		offspring[j] = psi[j] * r[j] + (1 - psi[j])*s[j];
	}
	return offspring;
}

//	Selection step of GA
template<typename T>
std::vector<T> Solver<T, GAstruct<T>>::selection()
{
	// Generate r and s indices
	T xi = quantile(dist, distribution(generator));
	size_t r = static_cast<size_t>(std::round(npop * xi));
	xi = quantile(dist, distribution(generator));
	size_t s = static_cast<size_t>(std::round(npop * xi));
	std::vector<T> offspring = crossover(individuals[r], individuals[s]);
	return offspring;
}

// Mutation step of GA
template<typename T>
std::vector<T> Solver<T, GAstruct<T>>::mutation(const std::vector<T>& individual)
{
	std::vector<T> mutated = individual;
	for (auto j = 0; j < ndv; ++j)
	{
		T r = distribution(generator);
		if (pi < r)
		{
			std::normal_distribution<T> ndistribution(0, stdev[j]);
			T epsilon = ndistribution(generator);
			mutated[j] = mutated[j] + epsilon;
		}
	}
	return mutated;
}

template<typename T>
template<typename F, typename C>
void Solver<T, GAstruct<T>>::run_algo(F f, C c)
{
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	find_min_cost(f);
	for (iter = 0; iter < iter_max; ++iter)
	{
		if (tol > std::abs(fitness_cost - opt))
		{
			break;
		}
		if (individuals.size() < 3)
		{
			break;
		}
		std::sort(individuals.begin(), individuals.end(), comparator);
		std::vector<std::vector<T>> offsprings;
		for (auto i = 0; i < npop; ++i)
		{
			std::vector<T> offspring = selection();
			offsprings.push_back(offspring);
		}
		npop = individuals.size();
		individuals.erase(individuals.begin() + nkeep, individuals.begin() + npop);
		npop = individuals.size();
		nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
		for (auto& p : individuals)
		{
			p = mutation(p);
			while (!c(p))
			{
				p = mutation(p);
			}
		}
		// Standard Deviation is not constant in GA
		for (auto j = 0; j < ndv; ++j)
		{
			stdev[j] = stdev[j] + 0.02 * stdev[j];
		}
		fitness_cost = f(individuals[0]);
	}
}