#pragma once

#include <boost/math/distributions.hpp>
#include "ealgorithm_base.h"

//! Genetic Algorithms (GA) Class
template<typename T>
class Genetic_Algo : public Solver_base<T, Genetic_Algo<T> >
{
public:
	//! Constructor
	Genetic_Algo(const T& i_x_rate, const T& i_pi, const T& i_alpha, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max) : x_rate{ i_x_rate }, pi{ i_pi }, alpha{ i_alpha }, stdev_mut{ i_stdev },
		Solver_base<T, Genetic_Algo<T> > { i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }
	{
		assert(x_rate > 0 && x_rate <= 1);
		assert(pi > 0 && pi <= 1);
		assert(alpha > 0);
	}
	//! Type of the algorithm :: string
	const std::string type = "Genetic Algorithms";
	//! Runs the algorithm until stopping criteria
	template<typename F, typename C> void run_algo(F f, C c);
private:
	//! Natural Selection rate
	const T x_rate;
	//! Probability of mutating
	const T pi;
	//! Parameter alpha for Beta distribution
	const T alpha;
	//! Number of individuals to be kept
	size_t nkeep;
	//! Standard Deviation is not constant in GA
	std::vector<T> stdev_mut;
	//! Beta distribution
	boost::math::beta_distribution<T> dist;
	//! Crossover step of GA
	std::vector<T> crossover(std::vector<T> r, std::vector<T> s);
	//!	Selection step of GA
	std::vector<T> selection();
	//! Mutation step of GA
	std::vector<T> mutation(const std::vector<T>& individual);
};

template<typename T>
std::vector<T> Genetic_Algo<T>::crossover(std::vector<T> r, std::vector<T> s)
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

template<typename T>
std::vector<T> Genetic_Algo<T>::selection()
{
	//! Generate r and s indices
	T xi = quantile(dist, distribution(generator));
	size_t r = static_cast<size_t>(std::round(individuals.size() * xi));
	xi = quantile(dist, distribution(generator));
	size_t s = static_cast<size_t>(std::round(individuals.size() * xi));
	//! Produce offsrping using r and s indices by crossover
	std::vector<T> offspring = crossover(individuals[r], individuals[s]);
	return offspring;
}

template<typename T>
std::vector<T> Genetic_Algo<T>::mutation(const std::vector<T>& individual)
{
	std::vector<T> mutated = individual;
	for (auto j = 0; j < ndv; ++j)
	{
		T r = distribution(generator);
		if (pi < r)
		{
			std::normal_distribution<T> ndistribution(0, stdev_mut[j]);
			T epsilon = ndistribution(generator);
			mutated[j] = mutated[j] + epsilon;
		}
	}
	return mutated;
}

template<typename T>
template<typename F, typename C>
void Genetic_Algo<T>
::run_algo(F f, C c)
{
	boost::math::beta_distribution<T> i_dist(1, alpha);
	dist = i_dist;
	nkeep = static_cast<size_t>(std::ceil(npop * x_rate));
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	start = std::chrono::system_clock::now();
	for (iter = 0; iter < iter_max; ++iter)
	{
		std::sort(individuals.begin(), individuals.end(), comparator);
		std::vector<std::vector<T>> offsprings;
		for (auto i = 0; i < individuals.size(); ++i)
		{
			std::vector<T> offspring = selection();
			offsprings.push_back(offspring);
		}
		//! Size of the population is mutable
		individuals.erase(individuals.begin() + nkeep, individuals.begin() + individuals.size());
		nkeep = static_cast<size_t>(std::ceil(individuals.size() * x_rate));
		for (auto& p : individuals)
		{
			p = mutation(p);
			while (!c(p))
			{
				p = mutation(p);
			}
		}
		//! Standard Deviation is not constant in GA
		for (auto& p : stdev_mut)
		{
			p = p + 0.02 * p;
		}
		fitness_cost = f(individuals[0]);
		if (tol > std::abs(fitness_cost))
		{
			break;
		}
		if (individuals.size() < 3)
		{
			break;
		}
	}
	end = std::chrono::system_clock::now();
}