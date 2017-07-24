#pragma once

#include "dependencies.h"

template<typename T>
class EA_base
{
public:
	void init_epsilon();
	void create_individuals(const std::vector<T>& decision_variables);
	void set_solver(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max);
	T get_tol() const { return tol };
	size_t get_iter_max() const { return iter_max; };
	size_t get_npop() const { return npop; };
	size_t get_ndv() const { return ndv; };
	std::vector<T> get_stdev() const { return stdev; };
	std::vector<std::vector<T>> get_individuals() const { return individuals; };
	//void set_standard_deviation(std::vector<T> stdev) { this.stdev = stdev; };
	//void set_population_size(std::vector<T> npop) { this.npop = npop; };
private:
	// Tolerance
	T tol;
	// Number of maximum iterations
	size_t iter_max;
	// Standard deviation of the decision variables
	std::vector<T> stdev;
	// Size of the population
	size_t npop;
	// Number of decision variables
	size_t ndv;
	// Population
	std::vector<std::vector<T>> individuals;
};

template<typename T>
void EA_base<T>::set_solver(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max)
{
	this->tol = tol;
	this->iter_max = iter_max;
	this->npop = npop;
	this->ndv = decision_variables.size();
	this->stdev = stdev;
	if (individuals.size() > 0)
	{
		individuals.erase(individuals.begin(), individuals.end());
	}
	create_individuals(decision_variables);
	init_epsilon();
}

template<typename T>
void EA_base<T>::init_epsilon()
{
	std::random_device generator;
	std::vector<T> diff(npop);
	//std::vector<T> sum(ndv);
	//std::vector<T> mean(ndv);
	//std::vector<T> stdev(ndv);
	//std::vector<T> epsilon(ndv);
	T epsilon;
	/*
	for (auto& p : sum)
	{
	p = 0.0;
	}
	for (auto& p : individuals)
	{
	for (auto j = 0; j < ndv; ++j)
	{
	sum[j] = sum[j] + p[j];
	}
	}
	for (auto j = 0; j < ndv; ++j)
	{
	mean[j] = sum[j] / (npop);
	}
	for (auto j = 0; j < ndv; ++j)
	{
	for (auto i = 0; i < npop; ++i)
	{
	diff[i] = individuals[i][j] - mean[j];
	}
	T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
	stdev[j] = std::sqrt(sq_sum / npop);
	}
	*/
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			std::normal_distribution<> ndistribution(0, stdev[j]);
			epsilon = ndistribution(generator);
			individuals[i][j] = individuals[i][j] + epsilon;
		}
	}
}

template<typename T>
void EA_base<T>::create_individuals(const std::vector<T>& decision_variables)
{
	individuals.resize(npop, std::vector<T>(ndv));
	for (auto i = 0; i < npop; ++i)
	{
		individuals[i] = decision_variables;
	}
}