#pragma once

#include "ealgorithm_base.h"

//! Differntial Evolution Structure, used in the actual algorithm and for type deduction
template<typename T>
struct DE : EA_base<T>
{
public:
	//! Constructor
	DE(const T& i_cr, const T& i_f_param, const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev, 
		const size_t& i_npop, const T& i_tol, const size_t& i_iter_max)
		: cr{ i_cr }, f_param{ i_f_param }, EA_base { i_decision_variables, i_stdev, i_npop, i_tol, i_iter_max }
	{
		assert(cr > 0 && cr <= 1);
		assert(f_param > 0 && f_param <= 1);
	}
	//! Crossover Rate
	const T cr;
	//! Mutation Scale Fuctor
	const T f_param;
};

//! Differential Evolution Algorithm Class
template<typename T>
class Solver<T, DE<T>> : public Solver_base<T>
{
public:
	//! Constructor
	template<typename F, typename C> Solver(const DE<T>& de, F f, C c) : Solver_base<T>{ { de.decision_variables, de.stdev, de.npop, de.tol, de.iter_max }, f, c }, cr{ de.cr }, f_param{ de.f_param }
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
	//! Type of the algorithm
	const std::string type = "Differential Evolution";
	//! Runs the algorithm until stopping criteria
	template<typename F, typename C> void run_algo(F f, C c);
private:
	//! Crossover Rate
	T cr;
	//! Mutation Scale Fuctor
	T f_param;
	//! Indices of population
	std::vector<size_t> indices;
	//! Random number engine for the indices
	std::mt19937 engine{ generator() };
	//! Uniform size_t distribution of the indices
	std::uniform_int_distribution<size_t> ind_distribution;
	//! Method that constructs the donor vector
	std::vector<T> construct_donor();
	//! Method that constructs the trial vector
	std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
};

template<typename T>
std::vector<T> Solver<T, DE<T>>::construct_donor()
{
	std::vector<T> donor(ndv);
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
	return donor;
}

template<typename T>
std::vector<T> Solver<T, DE<T>>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
{
	std::vector<T> trial(ndv);
	std::vector<size_t> j_indices;
	for (auto j = 0; j < ndv; ++j)
	{
		j_indices.push_back(j);
	}
	std::uniform_int_distribution<size_t> j_ind_distribution(0, ndv - 1);
	for (auto j = 0; j < ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = j_indices[j_ind_distribution(engine)];
		if (epsilon <= cr || j == jrand)
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

template<typename T>
template<typename F, typename C>
void Solver<T, DE<T>>::run_algo(F f, C c)
{
	//! Differential Evolution starts here
	for (iter = 0; iter < iter_max; ++iter)
	{
		//! Stopping Criteria
		if (tol > std::abs(fitness_cost))
		{
			break;
		}
		for (auto& p : individuals)
		{
			//! Construct donor and trial vectors
			std::vector<T> donor = construct_donor();
			while (!c(donor))
			{
				donor = construct_donor();
			}
			std::vector<T> trial = construct_trial(p, donor);
			if (f(trial) <= f(p))
			{
				p = trial;
			}
		}
		//! Recalculate minimum cost individual of the population
		find_min_cost(f);
	}
}