#pragma once

#include <boost/math/distributions.hpp>
#include "ealgorithm_base.h"

//! Replacing or remove individuals strategies during mutation
enum class Strategy { keep_same, re_mutate, remove };

//! Genetic Algorithms Structure, used in the actual algorithm and for type deduction
template<typename T>
struct GA : EA_base<T>
{
public:
	//! Constructor
	GA(const T& i_x_rate, const T& i_pi, const T& i_alpha, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max, const Strategy& i_strategy)
		: x_rate{ i_x_rate }, pi{ i_pi }, alpha{ i_alpha }, EA_base{ i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }, strategy { i_strategy }
	{
		assert(x_rate > 0 && x_rate <= 1);
		assert(pi > 0 && pi <= 1);
	}
	//!Natural Selection rate
	const T x_rate;
	//! Probability of mutating
	const T pi;
	//! Parameter alpha for Beta distribution
	const T alpha;
	Strategy strategy;
};

//! Genetic Algorithms (GA) Class
template<typename T, typename F, typename C>
class Solver<GA<T>, F, C> : public Solver_base<T, F, C>
{
public:
	//! Constructor
	Solver(const GA<T>& i_ga, F f, C c) : 
		Solver_base<T, F, C>{ i_ga.decision_variables, i_ga.npop, i_ga.stdev, i_ga.tol, f, c }, ga{ i_ga }, npop{ i_ga.npop }, stdev{ i_ga.stdev }
	{
		bdistribution = boost::math::beta_distribution<T>::beta_distribution(1, ga.alpha);
	}
	//! Type of the algorithm
	const std::string type = "Genetic Algorithms";
	//! Runs the algorithm until stopping criteria
	void run_algo();
private:
	//! Genetic Algorithms structure used internally
	const GA<T> ga;
	//! Size of the population is mutable
	size_t npop;
	//! Standard deviation is mutable, so a copy is created
	std::vector<T> stdev;
	//! Beta distribution
	boost::math::beta_distribution<T> bdistribution;
	//! Crossover step of GA
	std::vector<T> crossover(std::vector<T> r, std::vector<T> s);
	//!	Selection step of GA
	std::vector<T> selection();
	//! Mutation step of GA
	std::vector<T> mutation(const std::vector<T>& individual);
	//! Returns number of individuals to be kept in each generation, thus 
	size_t nkeep();
};

template<typename T, typename F, typename C>
std::vector<T> Solver<GA<T>, F, C>::crossover(std::vector<T> r, std::vector<T> s)
{
	std::vector<T> offspring(ga.ndv);
	std::vector<T> psi(ga.ndv);
	for (auto j = 0; j < ga.ndv; ++j)
	{
		psi[j] = distribution(generator);
		offspring[j] = psi[j] * r[j] + (1 - psi[j]) * s[j];
	}
	return offspring;
}

template<typename T, typename F, typename C>
std::vector<T> Solver<GA<T>, F, C>::selection()
{
	//! Generate r and s indices
	T xi = quantile(bdistribution, distribution(generator));
	size_t r = static_cast<size_t>(std::round(static_cast<T>(individuals.size()) * xi));
	xi = quantile(bdistribution, distribution(generator));
	size_t s = static_cast<size_t>(std::round(static_cast<T>(individuals.size()) * xi));
	//! Produce offsrping using r and s indices by crossover
	std::vector<T> offspring = crossover(individuals[r], individuals[s]);
	return offspring;
}

template<typename T, typename F, typename C>
std::vector<T> Solver<GA<T>, F, C>::mutation(const std::vector<T>& individual)
{
	std::vector<T> mutated = individual;
	for (auto j = 0; j < ga.ndv; ++j)
	{
		T r = distribution(generator);
		if (ga.pi < r)
		{
			std::normal_distribution<T> ndistribution(0, stdev[j]);
			T epsilon = ndistribution(generator);
			mutated[j] = mutated[j] + epsilon;
		}
	}
	return mutated;
}

template<typename T, typename F, typename C>
size_t Solver<GA<T>, F, C>::nkeep()
{
	return static_cast<size_t>(std::ceil(npop * ga.x_rate));
}

template<typename T, typename F, typename C>
void Solver<GA<T>, F, C>::run_algo()
{
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	for (iter = 0; iter < ga.iter_max; ++iter)
	{
		std::sort(individuals.begin(), individuals.end(), comparator);
		individuals.erase(individuals.begin() + nkeep(), individuals.begin() + individuals.size());
		min_cost = individuals[0];
		if (ga.tol > std::abs(f(min_cost)))
		{
			solved_flag = true;
			break;
		}
		if (individuals.size() > 500)
		{
			individuals.erase(individuals.begin() + 500, individuals.begin() + individuals.size());
		}
		for (auto i = 0; i < npop; ++i)
		{
			std::vector<T> offspring = selection();
			individuals.push_back(offspring);
		}
		for (auto i = 1; i < individuals.size(); ++i)
		{
			std::vector<T> mutated = mutation(individuals[i]);
			if (ga.strategy == Strategy::remove)
			{
				if (!c(mutated))
				{
					if (i == individuals.size() - 1)
					{
						individuals.pop_back();
					}
					else
					{
						individuals.erase(individuals.begin() + i, individuals.begin() + i + 1);
					}
				}
				else
				{
					individuals[i] = mutated;
				}
			}
			if (ga.strategy == Strategy::keep_same)
			{
				if (!c(mutated))
				{
					
				}
				else
				{
					individuals[i] = mutated;
				}
			}
			if (ga.strategy == Strategy::re_mutate)
			{
				while (!c(mutated))
				{
					mutated = mutation(individuals[i]);
				}
				individuals[i] = mutated;
			}
		}
		//! Set the new population size which previous population size + natural selection rate * population size
		npop = individuals.size();
		//! Standard Deviation is not constant in GA
		for (auto& p : stdev)
		{
			p = p + 0.02 * p;
		}
	}
}