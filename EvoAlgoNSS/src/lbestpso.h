/** \file lbestpso.h
* \author Ioannis Anagnostopoulos
* \brief Classes and functions for the lbest (Local Best) Particle Swarm Optimisation with a ring topology
*/

#pragma once

#include <unordered_map>
#include <boost/math/constants/constants.hpp>
#include "ealgorithm_base.h"
#include <array>

namespace ea
{
	/** \struct PSOl
	*  \brief Local Best Particle Swarm Optimisation Structure, used in the actual algorithm and for type deduction
	*/
	template<typename T>
	struct PSOl : EA_base<T>
	{
	public:
		/** \fn PSOl(const T& i_c, const T& i_w, const std::vector<T>& i_vmax, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max,
		const bool& i_use_penalty_method = false, const Constraints_type& i_constraints_type = Constraints_type::none,
		const bool& i_print_to_output = true, const bool& i_print_to_file = true)
		\brief Constructor
		\param i_c c parameter for velocity update
		\param i_w Inertia parameter for velocity update
		\param i_vmax Maximum velocity
		\param i_decision_variables The starting values of the decision variables
		\param i_stdev The standard deviation
		\param i_npop The population size
		\param i_tol The tolerance
		\param i_iter_max The maximum number of iterations
		\param i_use_penalty_method Whether to used penalties or not
		\param i_constraints_type What kind of constraints to use
		\param i_print_to_output Whether to print to terminal or not
		\param i_print_to_file Whether to print to a file or not
		\return A PSOl<T> object
		*/
		PSOl(const T& i_c, const T& i_w, const std::vector<T>& i_vmax, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
			const size_t& i_npop, const T& i_tol, const size_t& i_iter_max,
			const bool& i_use_penalty_method, const Constraints_type& i_constraints_type,
			const bool& i_print_to_output, const bool& i_print_to_file) :
			EA_base<T>(i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max, i_use_penalty_method, i_constraints_type, i_print_to_output, i_print_to_file),
			c(i_c),
			w(i_w),
			vmax(i_vmax)
		{
			assert(c > 0);
			assert(w > 0);
			for (const auto& p : vmax) { assert(p > 0); };
			assert(vmax.size() == this->ndv);
		}
		/** \brief Parameter c for velocity update */
		const T c;
		/** \brief Inertia Variant of PSO : Inertia */
		const T w;
		/** \brief  Velocity Clamping Variant of PSO : Maximum Velocity */
		const std::vector<T> vmax;
		/** \brief Type of the algorithm */
		const std::string type = "Local Best Particle Swarm Optimisation";
	};

	/*! \class Solver<PSOl, T, F, C>
	*  \brief Local Best Particle Swarm Optimisation (PSO) Class
	*/
	template<typename T, typename F, typename C>
	class Solver<PSOl, T, F, C> : public Solver_base<Solver<PSOl, T, F, C>, PSOl, T, F, C>
	{
	public:
		friend class Solver_base<Solver<PSOl, T, F, C>, PSOl, T, F, C>;
		/*! \fn Solver(const PSO<T>& i_pso, F f, C c)
		*  \brief Constructor
		*  \param i_pso The particle swarm optimisation parameter structure that is used to construct the solver
		*  \param f A reference to the objective function
		*  \param c A reference to the constraints function
		*  \return A Solver<PSO, T, F, C> object
		*/
		Solver(const PSOl<T>& i_pso, F f, C c) :
			Solver_base<Solver<PSOl, T, F, C>, PSOl, T, F, C>(i_pso, f, c),
			pso(this->solver_struct),
			w(i_pso.w),
			vmax(i_pso.vmax),
			nneigh(i_pso.npop),
			neighbours(set_neighbourhoods())
		{
			velocity.resize(pso.npop, std::vector<T>(pso.ndv));
			local_best.resize(nneigh);
			for (auto& p : velocity)
			{
				for (auto& n : p)
				{
					n = 0.0;
				}
			}
			for (const auto& p : this->individuals)
			{
				personal_best.push_back(p);
				personal_best_cost.push_back(f(p));
			}
			for (auto& p : local_best)
			{
				p = personal_best[0];
			}
			for (size_t i = 0; i < pso.npop; ++i)
			{
				for (const auto& index : neighbours[i])
				{
					if (personal_best_cost[index] < this->f(local_best[i]))
					{
						local_best[i] = personal_best[index];
					}
				}
			}
			find_min_local_best();
		}
	private:
		/** \brief Particle Swarm Optimisation structure used internally (reference to solver_struct) */
		const PSOl<T>& pso;
		/** \brief Inertia is mutable, so a copy is created */
		T w;
		/** \brief Maximum Velocity is mutable, so a copy is created */
		std::vector<T> vmax;
		/** \brief Personal best vector of the particles, holds the best position recorded for each particle */
		std::vector<std::vector<T>> personal_best;
		/** \brief Personal best cost vector of the particles */
		std::vector<T> personal_best_cost;
		/** \brief Local best vector, holds the best position recorded for each neighbourhood */
		std::vector<std::vector<T>> local_best;
		/** \brief Velocity of the particles */
		std::vector<std::vector<T>> velocity;
		/** \brief Number of neighbourhoods */
		const size_t nneigh;
		/** \brief Neighbours of each particle */
		std::unordered_map<size_t, std::array<size_t, 3>> neighbours;
		/** \fn set_neighbourhoods
		*  \brief Set the neighbourhoods of the algorithm using particle indices
		*  \return A map matching particle indices to neighbourhoods
		*/
		std::unordered_map<size_t, std::array<size_t, 3>> set_neighbourhoods();
		/** \fn position_update()
		*  \brief Position update of the particles
		*  \return void
		*/
		void position_update();
		/** \fn best_update()
		*  \brief This method sets the personal and local best solutions
		*  \return void
		*/
		void best_update();
		/** \fn find_min_local_best()
		*  \brief This is a faster way to calculate the minimum cost unless there is only one neighbourhood, in which case it is the same as find_min_cost(F f)
		*  \return void
		*/
		void find_min_local_best();
		/** \fn check_pso_criteria
		*  \brief Define the maximum radius stopping criterion
		*  \return true if criteria are met, false otherwise
		*/
		bool check_pso_criteria();
		/** \fn run_algo
		*  \brief Runs the algorithm until stopping criteria
		*  return void
		*/
		void run_algo();
		/** \fn euclid_distance
		*  \brief Euclidean Distance of two vectors
		*..\param x,y The two vectors for which the distance is calculated
		*  \return Distance as a floating-point number
		*/
		T euclid_distance(const std::vector<T>& x, const std::vector<T>& y)
		{
			T sum = 0;
			for (size_t i = 0; i < x.size(); ++i)
			{
				sum = sum + std::pow(x[i] - y[i], 2);
			}
			return std::sqrt(sum);
		}
		/*! \fn display_parameters()
		*  \brief Display PSO parameters
		*  \return A std::stringstream of the parameters
		*/
		std::stringstream display_parameters()
		{
			std::stringstream parameters;
			parameters << "C:" << "," << pso.c << ",";
			parameters << "Inertia:" << "," << pso.w << ",";
			parameters << "Maximum Velocity:" << "," << pso.vmax;
			return parameters;
		}
	};

	template<typename T, typename F, typename C>
	std::unordered_map<size_t, std::array<size_t, 3>> Solver<PSOl, T, F, C>::set_neighbourhoods()
	{
		std::unordered_map<size_t, std::array<size_t, 3>> neighbours;
		for (size_t i = 0; i < pso.npop; ++i)
		{
			if (i == 0)
			{
				neighbours[i] = { pso.npop - 1, i, i + 1 };
			}
			else
			{
				if (i == pso.npop - 1)
				{
					neighbours[i] = { pso.npop - 2, i, 0 };
				}
				else
				{
					neighbours[i] = { i - 1, i, i + 1 };
				}
			}
		}
		return neighbours;
	}

	template<typename T, typename F, typename C>
	void Solver<PSOl, T, F, C>::position_update()
	{
		for (size_t i = 0; i < pso.npop; ++i)
		{
			for (size_t j = 0; j < pso.ndv; ++j)
			{
				std::uniform_real_distribution<double> c(0, pso.c);
				velocity[i][j] = w * velocity[i][j] + c(generator) * (personal_best[i][j] - this->individuals[i][j])
					+ c(generator) * (local_best[i][j] - this->individuals[i][j]);
				if (velocity[i][j] > vmax[j])
				{
					velocity[i][j] = vmax[j];
				}
				this->individuals[i][j] = this->individuals[i][j] + velocity[i][j];
			}
			//! Checks that the candidate is feasible
			if (!this->c(this->individuals[i]))
			{
				this->individuals[i] = personal_best[i];
			}
			if (this->f(this->individuals[i]) < this->f(personal_best[i]))
			{
				personal_best[i] = this->individuals[i];
				personal_best_cost[i] = this->f(this->individuals[i]);
			}
			for (const auto& index : neighbours[i])
			{
				if (personal_best_cost[index] < this->f(local_best[i]))
				{
					local_best[i] = personal_best[index];
				}
			}
		}
	}

	template<typename T, typename F, typename C>
	void Solver<PSOl, T, F, C>::best_update()
	{
		for (size_t i = 0; i < pso.npop; ++i)
		{
		}
	}

	template<typename T, typename F, typename C>
	void Solver<PSOl, T, F, C>::find_min_local_best()
	{
		for (size_t k = 0; k < nneigh; ++k)
		{
			if (this->f(local_best[k]) < this->f(this->min_cost))
			{
				this->min_cost = local_best[k];
			}
		}
	}

	template<typename T, typename F, typename C>
	bool Solver<PSOl, T, F, C>::check_pso_criteria()
	{
		std::vector<T> distance(pso.npop);
		for (size_t i = 0; i < pso.npop; ++i)
		{
			distance[i] = euclid_distance(this->individuals[i], this->min_cost);
		}
		T rmax = distance[0];
		for (size_t i = 0; i < pso.npop; ++i)
		{
			if (rmax < distance[i])
			{
				rmax = distance[i];
			}
			else
			{
			}
		}
		if (pso.tol > std::abs(this->f(this->min_cost)))// || rmax < pso.tol)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	/** \brief Inverse square of pi constant */
	template<typename T>
	const double inv_pi_sq = 1 / std::pow(boost::math::constants::pi<T>(), 2);

	template<typename T, typename F, typename C>
	void Solver<PSOl, T, F, C>::run_algo()
	{
		//! Local Best Particle Swarm starts here
		for (size_t iter = 0; iter < pso.iter_max; ++iter)
		{
			position_update();
			//best_update();
			find_min_local_best();
			//! Inertia weight is updated - Linear
			//w = pso.w - (pso.w - 0.4) * (static_cast<T>(iter) / static_cast<T>(pso.iter_max));
			//! Non-linear
			w = pso.w - (pso.w - 0.4) * std::pow((static_cast<T>(iter) / static_cast<T>(pso.iter_max)), inv_pi_sq<T>);
			//w = 0.729;
			//w = 0.5 + distribution(generator) / 2;
			this->last_iter = iter;
			if (check_pso_criteria())
			{
				this->solved_flag = true;
				break;
			}
		}
	}

}