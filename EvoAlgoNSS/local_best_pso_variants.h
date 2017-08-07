#pragma once

#include "local_best_pso.h"

template<typename T>
struct PSOstruct_inertia final : PSOstruct<T>
{
	PSOstruct_inertia(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const T& i_w, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev,
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: PSOstruct{ i_c1 , i_c2 , i_sneigh, i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }, w{ i_w }
	{
		assert(w > 0);
	}
	const T w;
};

template<typename T>
struct PSOstruct_clamping final : PSOstruct<T>
{
	PSOstruct_clamping(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const T& i_alpha, const std::vector<T>& i_vmax, 
		const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: PSOstruct{ i_c1 , i_c2 , i_sneigh, i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }, alpha{ i_alpha }, vmax{ i_vmax }
	{
		assert(alpha > 0);
		for (const auto& p : vmax) { assert(p > 0); };
	}
	const T alpha;
	const std::vector<T> vmax;
};

// Inertia Variant of PSO
template<typename T>
class Solver<T, PSOstruct_inertia<T>> final : public Solver<T, PSOstruct<T>>
{
public:
	Solver(const PSOstruct_inertia<T>& pso) : Solver < T, PSOstruct<T>>{ { pso.c1, pso.c2, pso.sneigh, pso.decision_variables, pso.stdev, pso.npop, pso.tol, pso.iter_max} }, w{ pso.w }
	{
	}
	const std::string type = "Local Best Particle Swarm Optimisation with Inertia";
	template<typename F, typename C> void run_algo(F f, C c);
private:
	T w;
};

// Runs the algorithm until the stopping criteria
template<typename T>
template<typename F, typename C>
void Solver<T, PSOstruct_inertia<T>>::run_algo(F f, C c)
{
	find_min_cost(f);
	// Local Best Particle Swarm starts here
	for (iter = 0; iter < iter_max; ++iter)
	{
		check_pso_criteria();
		if (tol > std::abs(fitness_cost - opt) || rmax < tol)
		{
			break;
		}
		// Velocities of the particles is updated using inertia as well
		const auto& r = generate_r();
		for (auto i = 0; i < npop; ++i)
		{
			for (auto j = 0; j < ndv; ++j)
			{
				velocity[i][j] = w * velocity[i][j] + c1 * r[0][j] * (personal_best[i][j] - individuals[i][j])
					+ c2 * r[1][j] * (local_best[neighbourhoods[i]][j] - individuals[i][j]);
			}
		}
		position_update();
		check_particle_constraints(c);
		set_best(f);
		min_cost = find_min_local_best(f);
		w = ((w - 0.4) * (iter_max - iter)) / (iter_max + 0.4);
	}
}

// Velocity Clamping Variant of 
template<typename T>
class Solver<T, PSOstruct_clamping<T>> final : public Solver<T, PSOstruct<T>>
{
public:
	Solver(const PSOstruct_clamping<T>& pso) : Solver < T, PSOstruct<T>>{ { pso.c1, pso.c2, pso.sneigh, pso.decision_variables, pso.stdev, pso.npop, pso.tol, pso.iter_max} }
		, alpha{pso.alpha}, vmax{pso.vmax}
	{
		assert(vmax.size() == ndv);
	}
	const std::string type = "Local Best Particle Swarm Optimisation with Velocity Clamping";
	template<typename F, typename C> void run_algo(F f, C c);
private:
	T alpha;
	// Maximum Velocity
	std::vector<T> vmax;
};


// Runs the algorithm until the stopping criteria
template<typename T>
template<typename F, typename C>
void Solver<T, PSOstruct_clamping<T>>::run_algo(F f, C c)
{
	find_min_cost(f);
	// Local Best Particle Swarm starts here
	for (iter = 0; iter < iter_max; ++iter)
	{
		check_pso_criteria();
		if (tol > std::abs(fitness_cost - opt) || rmax < tol)
		{
			break;
		}
		velocity_update();
		// Velocities of the particles is updated using velocity clapming as well
		const auto& r = generate_r();
		for (auto i = 0; i < npop; ++i)
		{
			for (auto j = 0; j < ndv; ++j)
			{
				velocity[i][j] = velocity[i][j] + c1 * r[0][j] * (personal_best[i][j] - individuals[i][j])
					+ c2 * r[1][j] * (local_best[neighbourhoods[i]][j] - individuals[i][j]);
				if (velocity[i][j] > vmax[j])
				{
					velocity[i][j] = vmax[j];
				}
			}
		}
		position_update();
		check_particle_constraints(c);
		set_best(f);
		min_cost = find_min_local_best(f);
		// Maximum velocity is reduced using the current iteration and 
		for (auto& p : vmax)
		{
			p = (1 - std::pow(iter / iter_max, alpha)) * p;
		}
	}
}