#pragma once

#include "ealgorithm_base.h"

template<typename T>
struct DEstruct : EAstruct<T>
{
public:
	DEstruct(const T& i_cr, const T& i_f_param, const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: cr{ i_cr }, f_param{ i_f_param }, EAstruct{ i_npop, i_tol, i_iter_max }
	{
		assert(cr > 0 && cr <= 1);
		assert(f_param > 0 && f_param <= 1);
	}
	// Crossover Rate
	const T cr;
	// Mutation Scale Fuctor
	const T f_param;
};

//Differential Evolution Algorithm Class
template<typename T, typename F>
class Solver<T, F, DEstruct<T>> : public Solver<T, F, EAstruct<T>>
{
public:
	// Constructor for the Differential Evolution Class
	Solver(const DEstruct<T>& de, const Population<T>& popul) : Solver < T, F, EAstruct<T>>{ { de.npop, de.tol, de.iter_max }, popul }, cr{ de.cr }, f_param{ de.f_param }
	{
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		for (auto i = 0; i < npop; ++i)
		{
			indices.push_back(i);
		}
		std::uniform_int_distribution<size_t> i_ind_distribution(0, indices.size() - 1);
		ind_distribution = i_ind_distribution;
	}
private:
	// Crossover Rate
	T cr;
	// Mutation Scale Fuctor
	T f_param;
	// Indices of population
	std::vector<size_t> indices;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	std::random_device generator;
	// Uniform Real Distribution used for generating the random vectors for mutation and epsilon for crossover 
	std::uniform_real_distribution<T> distribution;
	std::uniform_int_distribution<size_t> ind_distribution;
	std::string get_type_of_solver() { return "Differential Evolution"; };
	void construct_donor(std::vector<T>& donor);
	void construct_trial(const std::vector<T>& target, const std::vector<T>& donor, std::vector<T>& trial);
	void modify_individuals(F f);
	void run_algo(F f, const T& opt);
};

// Method that constructs the donor vector
template<typename T, typename F>
void Solver<T, F, DEstruct<T>>::construct_donor(std::vector<T>& donor)
{
	std::vector<size_t> r_i;
	while (r_i.size() < 3)
	{
		r_i.push_back(indices[ind_distribution(engine)]);
		if (r_i.size() > 1 && r_i.end()[-1] == r_i.end()[-2])
		{
			r_i.pop_back();
		}
	}
	for (auto j = 0; j < ndv; ++j)
	{
		donor[j] = individuals[r_i[0]][j] + f_param * (individuals[r_i[1]][j] - individuals[r_i[2]][j]);
	}
}

// Method that constructs the trial vector
template<typename T, typename F>
void Solver<T, F, DEstruct<T>>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor, std::vector<T>& trial)
{
	for (auto j = 0; j < ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = indices[ind_distribution(engine)];
		if (epsilon <= cr || j == jrand)
		{
			trial[j] = donor[j];
		}
		else
		{
			trial[j] = target[j];
		}
	}
}

template<typename T, typename F>
void Solver<T, F, DEstruct<T>>::modify_individuals(F f)
{
	for (auto i = 0; i < npop; ++i)
	{
		std::vector<T> donor(ndv);
		std::vector<T> trial(ndv);
		// Construct donor and trial vectors
		construct_donor(donor);
		construct_trial(individuals[i], donor, trial);
		// Replace individual i with trial if trial's cost is lower / Selection step
		if (f(trial) <= f(individuals[i]))
		{
			individuals[i] = trial;
		}
	}
}

// Solve function overload for the Differential Evolution class
template<typename T, typename F>
void Solver<T,F, DEstruct<T>>::run_algo(F f, const T& opt)
{
	// Differential Evolution starts here
	for (auto iter = 0; iter < iter_max; ++iter)
	{
		// Stopping Criteria
		if (tol > std::abs(fitness_cost - opt))
		{
			last_iter = iter;
			break;
		}
		modify_individuals(f);
		// Recalculate minimum cost individual of the population
		find_min_cost(f);
		last_iter = iter;
	}
}