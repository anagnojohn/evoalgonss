#pragma once

#include "ealgorithm_base.h"
#include <type_traits>

// Euclidean Distance of two vectors
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

template<typename T>
struct PSOstruct : EAstruct<T>
{
public:
	PSOstruct(const T& i_c1, const T& i_c2, const size_t& i_sneigh, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, 
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: c1{ i_c1 }, c2{ i_c2 }, sneigh{ i_sneigh }, EAstruct{ i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }
	{
		assert(c1 > 0);
		assert(c2 > 0);
		assert(sneigh > 0);
		assert(sneigh < npop);
	}
	const T c1;
	const T c2;
	// Size of each neighbourhood
	const size_t sneigh;
};

// Local Best Particle Swarm Optimisation Class
template<typename T>
class Solver<T, PSOstruct<T>> : public Solver<T, EAstruct<T>>
{
public:
	Solver(const PSOstruct<T>& pso) : Solver <T, EAstruct<T>>{ { pso.decision_variables, pso.stdev, pso.npop, pso.tol, pso.iter_max } }, c1{ pso.c1 }, c2{ pso.c2 }, sneigh{ pso.sneigh }
	{
		init_pso();
	}
	const std::string type = "Local Best Particle Swarm Optimisation";
	template<typename F, typename C> void run_algo(F f, C c);
protected:
	T c1;
	T c2;
	// Neighbourhood size
	size_t sneigh;
	std::vector<std::vector<T>> personal_best;
	std::vector<std::vector<T>> local_best;
	std::vector<std::vector<T>> velocity;
	// Number of neighbourhoods
	size_t nneigh;
	// Neighbourhoods
	std::map<size_t, size_t> neighbourhoods;
	// Random number generator
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
	// Maximum Radius
	T rmax;
	// Distance betwwen individuals
	std::vector<T> distance;
	void init_pso();
	void velocity_update();
	void set_neighbourhoods();
	template<typename F> void set_best(F f);
	std::vector<std::vector<T>> generate_r();
	template<typename F> std::vector<T> find_min_local_best(F f);
	void position_update();
	void check_pso_criteria();
	template<typename C> void check_particle_constraints(C c);
};

// Initialisation of the algorithm
template<typename T>
void Solver<T, PSOstruct<T>>::init_pso()
{
	std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
	distribution = i_distribution;
	nneigh = static_cast<size_t>(std::ceil(npop / sneigh));
	velocity.resize(npop, std::vector<T>(ndv));
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
	for (auto k = 0; k < nneigh; ++k)
	{
		local_best.push_back(individuals[0]);
		local_best[k] = personal_best[0];
	}
	set_neighbourhoods();
	distance.resize(npop);
}

// Set the neighbourhoods of the algorithm
template<typename T>
void Solver<T, PSOstruct<T>>::set_neighbourhoods()
{
	size_t neigh_index = 0;
	size_t counter = 0;
	for (auto i = 0; i < npop; ++i)
	{
		
		if (counter < sneigh)
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

// This method generates r1 and r2 for the velocity update rule
template<typename T>
std::vector<std::vector<T>> Solver<T, PSOstruct<T>>::generate_r()
{
	std::vector<std::vector<T>> r(2, std::vector<T>(ndv));
	for (auto i = 0; i < 2; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			r[i][j] = (distribution(generator));
		}
	}
	return r;
}

// Velocity update of the particles
template<typename T>
void Solver<T, PSOstruct<T>>::velocity_update()
{
	const auto& r = generate_r();
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			velocity[i][j] = velocity[i][j] + c1 * r[0][j] * (personal_best[i][j] - individuals[i][j])
				+ c2 * r[1][j] * (local_best[neighbourhoods[i]][j] - individuals[i][j]);
		}
	}
}

// Position update of the particles
template<typename T>
void Solver<T, PSOstruct<T>>::position_update()
{
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			{
				individuals[i][j] = individuals[i][j] + velocity[i][j];	
			}
		}
	}
}

// This method sets the personal and local best solutions
template<typename T>
template<typename F>
void Solver<T, PSOstruct<T>>::set_best(F f)
{
	for (auto i = 0; i < npop; ++i)
	{
		if (f(individuals[i]) < f(personal_best[i]))
		{
			personal_best[i] = individuals[i];
		}
	}
	for (auto i = 0; i < npop; ++i)
	{
		if (f(personal_best[i]) < f(local_best[neighbourhoods[i]]))
		{
			local_best[neighbourhoods[i]] = personal_best[i];
		}
	}
}

// This is a faster way to calculate the minimum cost unless there is only one neighbourhood, in which case it is the same as find_min_cost(F f)
template<typename T>
template<typename F>
std::vector<T> Solver<T, PSOstruct<T>>::find_min_local_best(F f)
{
	std::vector<T> min_cost(ndv);
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

// Define the maximum radius stopping criterion
template<typename T>
void Solver<T, PSOstruct<T>>::check_pso_criteria()
{
	for (auto i = 0; i < npop; ++i)
	{
		distance[i] = euclid_distance(individuals[i], min_cost);
	}
	rmax = distance[0];
	for (auto i = 0; i < npop; ++i)
	{
		if (rmax < distance[i])
		{
			rmax = distance[i];
		}
	}
}

template<typename T>
template<typename C>
void Solver<T, PSOstruct<T>>::check_particle_constraints(C c)
{
	for (auto i = 0; i < npop; ++i)
	{
		if (!c(individuals[i]))
		{
			individuals[i] = personal_best[i];
		}
	}
}

// Runs the algorithm until the stopping criteria
template<typename T>
template<typename F, typename C>
void Solver<T, PSOstruct<T>>::run_algo(F f, C c)
{
	for (auto i = 0; i < npop; ++i)
	{
		if (f(personal_best[i]) < f(local_best[neighbourhoods[i]]))
		{
			local_best[neighbourhoods[i]] = personal_best[i];
		}
	}
	find_min_local_best(F f);
	// Local Best Particle Swarm starts here
	for (iter = 0; iter < iter_max; ++iter)
	{	
		check_pso_criteria();
		if (tol > std::abs(fitness_cost - opt) || rmax < tol)
		{
			break;
		}
		velocity_update();
		position_update();
		check_particle_constraints(c);
		set_best(f);
		min_cost = find_min_local_best(f);
	}
}