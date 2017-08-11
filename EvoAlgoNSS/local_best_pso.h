#pragma once

#include "ealgorithm_base.h"
#include <unordered_map>

//! Particle Swarm Optimisation Structure, used in the actual algorithm and for type deduction
template<typename T>
struct PSO : EA_base<T>
{
public:
	//! Constructor
	PSO(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const T& i_w, const T& i_alpha, const std::vector<T>& i_vmax, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: c1{ i_c1 }, c2{ i_c2 }, sneigh{ i_sneigh }, EA_base{ i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }, w{ i_w }, alpha{ i_alpha }, vmax{ i_vmax }
	{
		assert(c1 > 0);
		assert(c2 > 0);
		assert(sneigh > 0);
		assert(sneigh < npop);
		assert(w > 0);
		assert(alpha > 0);
		for (const auto& p : vmax) { assert(p > 0); };
		assert(vmax.size() == ndv);
	}
	//! Parameter c1 for velocity update
	const T c1;
	//! Parameter c2 for velocity update
	const T c2;
	//! Neighbourhood size
	const size_t sneigh;
	//! Inertia Variant of PSO : Inertia
	const T w;
	//! Alpha Parameter for maximum velocity
	const T alpha;
	//! Velocity Clamping Variant of PSO : Maximum Velocity
	const std::vector<T> vmax;
};

//! Local Best Particle Swarm Optimisation Class
template<typename T>
class Solver<PSO<T>> : public Solver_base<T>
{
public:
	//! Constructor
	template<typename F, typename C> Solver(const PSO<T>& i_pso, F f, C c) : Solver_base<T>{ i_pso.decision_variables, i_pso.npop, i_pso.stdev, f, c }, w{ i_pso.w }, vmax{ i_pso.vmax }, pso{ i_pso }
	{
		nneigh = static_cast<size_t>(std::ceil(pso.npop / pso.sneigh));
		velocity.resize(pso.npop, std::vector<T>(pso.ndv));
		local_best.resize(nneigh);
		for (auto& p : velocity)
		{
			for (auto& n : p)
			{
				n = 0.0;
			}
		}
		for (const auto& p : individuals)
		{
			personal_best.push_back(p);
		}
		for (auto& p : local_best)
		{
			p = personal_best[0];
		}
		set_neighbourhoods();
		distance.resize(pso.npop);
		for (auto i = 0; i < pso.npop; ++i)
		{
			if (f(personal_best[i]) < f(local_best[neighbourhoods[i]]))
			{
				local_best[neighbourhoods[i]] = personal_best[i];
			}
		}
		find_min_local_best(f);
	}
	//! Type of the algorithm
	const std::string type = "Local Best Particle Swarm Optimisation";
	//! Runs the algorithm until stopping criteria
	template<typename F, typename C> void run_algo(F f, C c);
protected:
	//! Particle Swarm Optimisation structure used internally
	const PSO<T>& pso;
	//! Inertia is mutable, so a copy is created
	T w;
	//! Maximum Velocity is mutable, so a copy is created
	std::vector<T> vmax;
	//! Personal best vector of the particles, holds the best position recorded for each particle
	std::vector<std::vector<T>> personal_best;
	//! Local best vector, holds the best position recorded for each neighbourhood
	std::vector<std::vector<T>> local_best;
	//! Velocity of the particles
	std::vector<std::vector<T>> velocity;
	//! Number of neighbourhoods
	size_t nneigh;
	//! Neighbourhoods
	std::unordered_map<size_t, size_t> neighbourhoods;
	//! Maximum Radius
	T rmax;
	//! Distance betwwen individuals
	std::vector<T> distance;
	//! Set the neighbourhoods of the algorithm using particle indices
	void set_neighbourhoods();
	//! This method generates r1 and r2 for the velocity update rule
	std::vector<std::vector<T>> generate_r();
	//! Velocity update of the particles
	void velocity_update();
	//! Position update of the particles
	void position_update();
	//! This method sets the personal and local best solutions
	template<typename F> void best_update(F f);
	//! This is a faster way to calculate the minimum cost unless there is only one neighbourhood, in which case it is the same as find_min_cost(F f)
	template<typename F> std::vector<T> find_min_local_best(F f);
	//! Define the maximum radius stopping criterion
	void check_pso_criteria();
	//! Checks that the candidates are feasible
	template<typename C> void check_particle_constraints(C c);
	//! Euclidean Distance of two vectors
	template<typename T>
	T euclid_distance(const std::vector<T>& x, const std::vector<T>& y)
	{
		T sum = 0;
		for (auto i = 0; i < x.size(); ++i)
		{
			sum = sum + std::pow(x[i] - y[i], 2);
		}
		return std::sqrt(sum);
	}
};

template<typename T>
void Solver<PSO<T>>::set_neighbourhoods()
{
	size_t neigh_index = 0;
	size_t counter = 0;
	for (auto i = 0; i < pso.npop; ++i)
	{
		if (counter < pso.sneigh)
		{
			neighbourhoods[i] = neigh_index;
			counter = counter + 1;
		}
		else
		{
			neigh_index = neigh_index + 1;
			counter = 0;
		}
	}
}

template<typename T>
std::vector<std::vector<T>> Solver<PSO<T>>::generate_r()
{
	std::vector<std::vector<T>> r(2, std::vector<T>(pso.ndv));
	for (auto i = 0; i < 2; ++i)
	{
		for (auto j = 0; j < pso.ndv; ++j)
		{
			r[i][j] = (distribution(generator));
		}
	}
	return r;
}

template<typename T>
void Solver<PSO<T>>::velocity_update()
{
	const auto& r = generate_r();
	for (auto i = 0; i < pso.npop; ++i)
	{
		for (auto j = 0; j < pso.ndv; ++j)
		{
			velocity[i][j] = w * velocity[i][j] + pso.c1 * r[0][j] * (personal_best[i][j] - individuals[i][j])
				+ pso.c2 * r[1][j] * (local_best[neighbourhoods[i]][j] - individuals[i][j]);
			if (velocity[i][j] > vmax[j])
			{
				velocity[i][j] = vmax[j];
			}
		}
	}
}

template<typename T>
void Solver<PSO<T>>::position_update()
{
	for (auto i = 0; i < pso.npop; ++i)
	{
		for (auto j = 0; j < pso.ndv; ++j)
		{
			individuals[i][j] = individuals[i][j] + velocity[i][j];
		}
	}
}

template<typename T>
template<typename F>
void Solver<PSO<T>>::best_update(F f)
{
	for (auto i = 0; i < pso.npop; ++i)
	{
		if (f(individuals[i]) < f(personal_best[i]))
		{
			personal_best[i] = individuals[i];
		}
	}
	for (auto i = 0; i < pso.npop; ++i)
	{
		if (f(personal_best[i]) < f(local_best[neighbourhoods[i]]))
		{
			local_best[neighbourhoods[i]] = personal_best[i];
		}
	}
}

template<typename T>
template<typename F>
std::vector<T> Solver<PSO<T>>::find_min_local_best(F f)
{
	std::vector<T> min_cost(pso.ndv);
	min_cost = local_best[0];
	for (auto k = 0; k < nneigh; ++k)
	{
		if (f(local_best[k]) < f(min_cost))
		{
			min_cost = local_best[k];
		}
	}
	fitness_cost = f(min_cost);
	return min_cost;
}

template<typename T>
void Solver<PSO<T>>::check_pso_criteria()
{
	for (auto i = 0; i < pso.npop; ++i)
	{
		distance[i] = euclid_distance(individuals[i], min_cost);
	}
	rmax = distance[0];
	for (auto i = 0; i < pso.npop; ++i)
	{
		if (rmax < distance[i])
		{
			rmax = distance[i];
		}
		else
		{

		}
	}
}


template<typename T>
template<typename C>
void Solver<PSO<T>>::check_particle_constraints(C c)
{
	for (auto i = 0; i < pso.npop; ++i)
	{
		if (!c(individuals[i]))
		{
			individuals[i] = personal_best[i];
		}
		else
		{

		}
	}
}

template<typename T>
template<typename F, typename C>
void Solver<PSO<T>>::run_algo(F f, C c)
{
	//! Local Best Particle Swarm starts here
	for (iter = 0; iter < pso.iter_max; ++iter)
	{
		velocity_update();
		position_update();
		check_particle_constraints(c);
		best_update(f);
		min_cost = find_min_local_best(f);
		//! Velocities of the particles is updated using inertia as well
		w = ((w - 0.4) * (pso.iter_max - iter)) / (pso.iter_max + 0.4);
		//! Maximum velocity is reduced using the current iteration and 
		for (auto& p : vmax)
		{
			p = (1 - std::pow(iter / pso.iter_max, pso.alpha)) * p;
		}
		check_pso_criteria();
		if (pso.tol > std::abs(fitness_cost) || rmax < pso.tol)
		{
			break;
		}
	}
}