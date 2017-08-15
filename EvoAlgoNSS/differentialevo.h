#pragma once

#include "ealgorithm_base.h"

namespace ea
{
	//! Differntial Evolution Structure, used in the actual algorithm and for type deduction
	template<typename T>
	struct DE : EA_base<T>
	{
	public:
		//! Constructor
		DE(const T& i_cr, const T& i_f_param, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
			const size_t& i_npop, const T& i_tol, const size_t& i_iter_max) :
			EA_base<T>{ i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max },
			cr{ i_cr },
			f_param{ i_f_param }
		{
			assert(cr > 0 && cr <= 1);
			assert(f_param > 0 && f_param <= 1);
		}
		//! Crossover Rate
		const T cr;
		//! Mutation Scale Fuctor
		const T f_param;
	};

	//! Differential Evolution Algorithm (DE) Class
	template<typename T, typename F, typename C>
	class Solver<DE, T, F, C> : public Solver_base<Solver<DE, T, F, C>, DE, T, F, C>
	{
	public:
		//! Constructor
		Solver(const DE<T>& i_de, const F& f, const C& c) :
                Solver_base<Solver<DE, T, F, C>, DE, T, F, C>{ i_de, f, c },
                de{ this->solver_struct },
                indices{ set_indices() }
		{
            ind_distribution = std::uniform_int_distribution<size_t>(0,de.npop - 1);
        };
		//! Type of the algorithm :: string
		const std::string type = "Differential Evolution";
		//! Method that displays the parameters of DE
		std::stringstream display_parameters()
		{
			std::stringstream parameters;
			parameters << "Crossover Rate: " << de.cr << "\n";
			parameters << "Mutation Scale Factor: " << de.f_param << "\n";
			return parameters;
		}
		//! Runs the algorithm until stopping criteria
		void run_algo();
	private:
		//! Differential Evolution structure used internally (reference to solver_struct)
		const DE<T>& de;
		//! Indices of population
		const std::vector<size_t> indices;
		//! Generate the indices
		std::vector<size_t> set_indices()
		{
			std::vector<size_t> indices;
			for (size_t i = 0; i < de.npop; ++i)
			{
				indices.push_back(i);
			}
			return indices;
		};
		//! Uniform size_t distribution of the indices
		std::uniform_int_distribution<size_t> ind_distribution;
		//! Method that constructs the donor vector
		std::vector<T> construct_donor();
		//! Method that constructs the trial vector
		std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
	};

	template<typename T, typename F, typename C>
	std::vector<T> Solver<DE, T, F, C>::construct_donor()
	{
		std::vector<T> donor(de.ndv);
		std::vector<size_t> r_i;
		//! Check that the indices are not the same
		while (r_i.size() < 3)
		{
			r_i.push_back(indices[ind_distribution(generator)]);
			if (r_i.size() > 1 && r_i.end()[-1] == r_i.end()[-2])
			{
				r_i.pop_back();
			}
		}
		for (size_t j = 0; j < de.ndv; ++j)
		{
			donor[j] = this->individuals[r_i[0]][j] + de.f_param * (this->individuals[r_i[1]][j] - this->individuals[r_i[2]][j]);
		}
		return donor;
	}

	template<typename T, typename F, typename C>
	std::vector<T> Solver<DE, T, F, C>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
	{
		std::vector<T> trial(de.ndv);
		std::vector<size_t> j_indices(de.ndv);
		for (size_t j = 0; j < de.ndv; ++j)
		{
			j_indices[j] = j;
		}
		std::uniform_int_distribution<size_t> j_ind_distribution(0, de.ndv - 1);
		for (size_t j = 0; j < de.ndv; ++j)
		{
			const T& epsilon = this->distribution(generator);
			const size_t& jrand = j_indices[j_ind_distribution(generator)];
			if (epsilon <= de.cr || j == jrand)
			{
				trial[j] = donor[j];
			}
			else
			{
				trial[j] = target[j];
			}
		}
		return trial;
	}

	template<typename T, typename F, typename C>
	void Solver<DE, T, F, C>::run_algo()
	{
		//! Differential Evolution starts here
		for (size_t iter = 0; iter < de.iter_max; ++iter)
		{
			for (auto& p : this->individuals)
			{
				//! Construct donor and trial vectors
				std::vector<T> donor = construct_donor();
				while (!this->c(donor))
				{
					donor = construct_donor();
				}
				const std::vector<T>& trial = construct_trial(p, donor);
				if (this->f(trial) <= this->f(p))
				{
					p = trial;
				}
			}
			//! Recalculate minimum cost individual of the population
			this->find_min_cost();
			//! Stopping Criteria
			if (de.tol > std::abs(this->f(this->min_cost)))
			{
				this->solved_flag = true;
				this->last_iter = iter;
				break;
			}
			this->last_iter = de.iter_max;
		}
	}
}