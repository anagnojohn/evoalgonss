#pragma once

#include "ealgorithm_base.h"
#include <boost/math/distributions.hpp>

namespace ea
{
	//! Replacing or remove individuals strategies during mutation
	enum class Strategy { keep_same, re_mutate, remove, none };

	//! Genetic Algorithms Structure, used in the actual algorithm and for type deduction
	template<typename T>
	struct GA : EA_base<T>
	{
	public:
		//! Constructor
		GA(const T& i_x_rate, const T& i_pi, const T& i_alpha, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
			const size_t& i_npop, const T& i_tol, const size_t& i_iter_max, const bool& i_use_penalty_method = false,
			const Constraints_type& i_constraints_type = Constraints_type::none, const Strategy& i_strategy = Strategy::keep_same,
			const bool& i_print_to_output = true, const bool& i_print_to_file = true)
			:  EA_base<T> { i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max, i_use_penalty_method, i_constraints_type, i_print_to_output, i_print_to_file },
			   x_rate{ i_x_rate }, pi{ i_pi }, alpha{ i_alpha }, strategy{ i_strategy }
		{
			assert(x_rate > 0 && x_rate <= 1);
			assert(pi > 0 && pi <= 1);
		}
		//! Natural Selection rate
		const T x_rate;
		//! Probability of mutating
		const T pi;
		//! Parameter alpha for Beta distribution
		const T alpha;
		//! Replacing or remove individuals strategies during mutation
		const Strategy strategy;
		//! Type of the algorithm
		const std::string type = "Genetic Algorithms";
	};

	//! Genetic Algorithms (GA) Class
	template<typename T, typename F, typename C>
	class Solver<GA, T, F, C> : public Solver_base<Solver<GA, T, F, C>, GA, T, F, C>
	{
	public:
		friend class Solver_base<Solver<GA, T, F, C>, GA, T, F, C>;
		//! Constructor
		Solver(const GA<T>& i_ga, F f, C c) :
			Solver_base<Solver<GA, T, F, C>, GA, T, F, C> { i_ga, f, c }, ga{ this->solver_struct }, npop{ i_ga.npop }, stdev{ i_ga.stdev }
		{
			bdistribution = boost::math::beta_distribution<T>(1, ga.alpha);
		}
	private:
		//! Genetic Algorithms structure used internally (reference to solver_struct)
		const GA<T>& ga;
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
		//! Runs the algorithm until stopping criteria
		void run_algo();
		//! Display GA parameters
		std::stringstream display_parameters()
		{
			std::stringstream parameters;
			parameters << "Natural Selection Rate" << "," << ga.x_rate << ",";
			parameters << "Probability of Mutation" << "," << ga.pi << ",";
			parameters << "Beta Distribution alpha" << "," << ga.alpha << ",";
			parameters << "Strategy" << ",";
			switch (ga.strategy)
			{
			case Strategy::keep_same: parameters << "Keep same individual"; break;
			case Strategy::re_mutate: parameters << "Re-mutate individual"; break;
			case Strategy::remove: parameters << "Remove individual"; break;
			default: parameters << "Do nothing"; break;
			}
			return parameters;
		}
	};

	template<typename T, typename F, typename C>
	std::vector<T> Solver<GA, T, F, C>::crossover(std::vector<T> r, std::vector<T> s)
	{
		std::vector<T> offspring(ga.ndv);
		std::vector<T> psi(ga.ndv);
		for (size_t j = 0; j < ga.ndv; ++j)
		{
			psi[j] = this->distribution(generator);
			offspring[j] = psi[j] * r[j] + (1 - psi[j]) * s[j];
		}
		return offspring;
	}

	template<typename T, typename F, typename C>
	std::vector<T> Solver<GA, T, F, C>::selection()
	{
		//! Generate r and s indices
		T xi = quantile(bdistribution, this->distribution(generator));
		size_t r = static_cast<size_t>(std::round(static_cast<T>(this->individuals.size()) * xi));
		xi = quantile(bdistribution, this->distribution(generator));
		size_t s = static_cast<size_t>(std::round(static_cast<T>(this->individuals.size()) * xi));
		//! Produce offsrping using r and s indices by crossover
		std::vector<T> offspring = crossover(this->individuals[r], this->individuals[s]);
		return offspring;
	}

	template<typename T, typename F, typename C>
	std::vector<T> Solver<GA, T, F, C>::mutation(const std::vector<T>& individual)
	{
		std::vector<T> mutated = individual;
		for (size_t j = 0; j < ga.ndv; ++j)
		{
			const T r = this->distribution(generator);
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
	size_t Solver<GA, T, F, C>::nkeep()
	{
		return static_cast<size_t>(std::ceil(this->npop * ga.x_rate));
	}

	template<typename T, typename F, typename C>
	void Solver<GA, T, F, C>::run_algo()
	{
		auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
		{
			return this->f(l) < this->f(r);
		};
		for (size_t iter = 0; iter < ga.iter_max; ++iter)
		{
			std::sort(this->individuals.begin(), this->individuals.end(), comparator);
			this->individuals.erase(this->individuals.begin() + nkeep(), this->individuals.begin() + this->individuals.size());
			this->min_cost = this->individuals[0];
			if (ga.tol > std::abs(this->f(this->min_cost)))
			{
				this->solved_flag = true;
				this->last_iter = iter;
				break;
			}
			if (this->individuals.size() > 1000)
			{
				//this->individuals.erase(this->individuals.begin() + 1000, this->individuals.begin() + this->individuals.size());
			}
			npop = this->individuals.size();
			for (size_t i = 0; i < npop; ++i)
			{
				std::vector<T> offspring = selection();
				this->individuals.push_back(offspring);
			}
			for (size_t i = 1; i < this->individuals.size(); ++i)
			{
				std::vector<T> mutated = mutation(this->individuals[i]);
				if (!this->c(mutated))
				{
					switch (ga.strategy)
					{
					case Strategy::keep_same: mutated = this->individuals[i]; break;
					case Strategy::re_mutate: 
					{
						while (!this->c(mutated))
						{
							mutated = mutation(this->individuals[i]);
						}
						break;
					}
					case Strategy::remove:
					{
						if (i == this->individuals.size() - 1)
						{
							this->individuals.pop_back();
						}
						else
						{
							this->individuals.erase(this->individuals.begin() + i, this->individuals.begin() + i + 1);
						}
						break;
					}
					case Strategy::none: break;
					}
				}
				this->individuals[i] = mutated;
			}
			//! Set the new population size which previous population size + natural selection rate * population size
			this->npop = this->individuals.size();
			//! Standard Deviation is not constant in GA
			for (auto& p : stdev)
			{
				p = p + 0.02 * p;
			}
		}
		this->last_iter = ga.iter_max;
	}
	
}