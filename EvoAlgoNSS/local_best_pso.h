#pragma once

#include "ealgorithm_base.h"

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
template<typename T, typename F>
class Solver<T, F, PSOstruct<T>> : public Solver<T,F, EAstruct<T>>
{
public:
	Solver(const PSOstruct<T>& pso) : Solver < T, F, EAstruct<T>>{ { pso.decision_variables, pso.stdev, pso.npop, pso.tol, pso.iter_max } }, c1{ pso.c1 }, c2{ pso.c2 }, sneigh{ pso.sneigh }
	{
		init_pso();
	}
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
	std::vector<std::vector<size_t>> neighbourhoods;
	// Indices of the population
	std::vector<size_t> indices;
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
	// Maximum Radius
	T rmax;
	// Distance betwwen individuals
	std::vector<T> distance;
	std::string get_type_of_solver() { return "Local Best Particle Swarm Optimisation"; };
	void init_pso();
	virtual void velocity_update(const size_t& iter);
	void set_neighbourhoods();
	void set_local_best(F f);
	std::vector<std::vector<T>> generate_r();
	std::vector<T> find_min_local_best(F f);
	void position_update();
	void run_algo(F f, const T& opt);
};

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::init_pso()
{
	std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
	distribution = i_distribution;
	nneigh = static_cast<size_t>(std::ceil(npop / sneigh));
	neighbourhoods.resize(nneigh);
	velocity.resize(npop, std::vector<T>(ndv));
	local_best.resize(nneigh);
	for (auto& p : velocity)
	{
		for (auto& n : p)
		{
			n = 0.0;
		}
	}
	for (auto k = 0; k < neighbourhoods.size(); ++k)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			local_best[k].push_back(individuals[0][j]);
		}
	}
	for (const auto& p : individuals)
	{
		personal_best.push_back(p);
	}
	for (auto i = 0; i < individuals.size(); ++i)
	{
		indices.push_back(i);
	}
	set_neighbourhoods();
	distance.resize(npop);
}

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::set_neighbourhoods()
{
	const auto& nneigh = neighbourhoods.size();
	size_t counter = 0;
	for (auto k = 0; k < nneigh; ++k)
	{
		while (counter < sneigh)
		{
			neighbourhoods[k].push_back(indices[k * sneigh + counter]);
			counter = counter + 1;
		}
	}
}

template<typename T, typename F>
std::vector<std::vector<T>> Solver<T, F, PSOstruct<T>>::generate_r()
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

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::velocity_update(const size_t& iter)
{
	const auto& r = generate_r();
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			for (auto j = 0; j < ndv; ++j)
			{
				{
					velocity[i][j] = velocity[i][j] + c1 * r[0][j] * (personal_best[i][j]
						- individuals[i][j]) + c2 * r[1][j] * (local_best[k][j] - individuals[i][j]);
				}
			}
		}
	}
}

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::position_update()
{
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			for (auto j = 0; j < ndv; ++j)
			{
				{
					individuals[i][j] = individuals[i][j] + velocity[i][j];
				}
			}
		}
	}
}

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::set_local_best(F f)
{
	for (auto i = 0; i < npop; ++i)
	{
		if (f(individuals[i]) < f(personal_best[i]))
		{
			personal_best[i] = individuals[i];
		}
	}
	for (auto k = 0; k < nneigh; ++k)
	{
		for (auto l = 0; l < neighbourhoods[k].size(); ++l)
		{
			auto i = neighbourhoods[k][l];
			if (f(personal_best[i]) < f(local_best[k]))
			{
				for (auto j = 0; j < ndv; j++)
				{
					local_best[k][j] = personal_best[i][j];
				}
			}
		}
	}
}


template<typename T, typename F>
std::vector<T> Solver<T, F, PSOstruct<T>>::find_min_local_best(F f)
{
	std::vector<T> min_cost(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		min_cost[j] = local_best[0][j];
	}
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

template<typename T, typename F>
void Solver<T, F, PSOstruct<T>>::run_algo(F f, const T& opt)
{
	find_min_cost(f);
	// Local Best Particle Swarm starts here
	for (auto iter = 0; iter < iter_max; ++iter)
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
		if (tol > std::abs(fitness_cost - opt) || rmax < tol)
		{
			last_iter = iter;
			break;
		}
		velocity_update(iter);
		position_update();
		set_local_best(f);
		min_cost = find_min_local_best(f);
		last_iter = iter;
	}
}