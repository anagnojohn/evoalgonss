#pragma once

#include <vector>

template<typename T>
class EAparams
{
public:
	EAparams(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max, const std::vector<T>& stdev)
	{
		set_params(decision_variables, npop, tol, iter_max, stdev);
	}
	// Getters for parameters
	T get_tol() const { return tol; };
	size_t get_iter_max() const { return iter_max; };
	size_t get_npop() const { return npop; ; };
	size_t get_ndv() const { return ndv; };
	std::vector<T> get_stdev() const { return stdev; };
	std::vector<std::vector<T>> get_individuals() const { return individuals; };
	// Setters for parameters
	void set_tol(const T& tol) { assert(tol > 0); this->tol = tol; };
	void set_npop(const size_t& npop) { assert(npop > 0); this->npop = npop; create_individuals(); };
	void set_iter_max(const size_t& iter_max) { assert(iter_max > 0); this->iter_max = iter_max; };
	void set_stdev(const std::vector<T>& stdev)
	{
		assert(decision_variables.size() == stdev.size());
		for (const auto& p : stdev)
		{
			assert(p > 0);
		}
		this->stdev = stdev;
		create_individuals;
		init_epsilon();
	};
	void set_decision_variables(const std::vector<T>& decision_variables)
	{
		this->decision_variables = decision_variables; ndv = decision_variables.size(); create_individuals();
	};
private:
	// Tolerance
	T tol;
	// Number of maximum iterations
	size_t iter_max;
	// Standard deviation of the decision variables
	std::vector<T> stdev;
	// Size of the population
	size_t npop;
	// Decision Variables
	std::vector<T> decision_variables;
	// Number of decision variables
	size_t ndv;
	// Population
	std::vector<std::vector<T>> individuals;
	void init_epsilon();
	void create_individuals();
	void set_params(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max, const std::vector<T>& stdev);
};

template<typename T>
void EAparams<T>::set_params(const std::vector<T>& decision_variables, const size_t& npop, const T& tol, const size_t& iter_max, const std::vector<T>& stdev)
{
	assert(decision_variables.size() > 0);
	assert(decision_variables.size() == stdev.size());
	assert(tol > 0);
	assert(iter_max > 0);
	assert(npop > 0);
	for (const auto& p : stdev)
	{
		assert(p > 0);
	}
	this->decision_variables = decision_variables;
	this->tol = tol;
	this->iter_max = iter_max;
	this->npop = npop;
	this->ndv = decision_variables.size();
	this->stdev = stdev;
	if (individuals.size() > 0)
	{
		individuals.erase(individuals.begin(), individuals.end());
	}
	create_individuals();
	init_epsilon();
}

template<typename T>
void EAparams<T>::init_epsilon()
{
	std::random_device generator;
	T epsilon;
	
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
void EAparams<T>::create_individuals()
{
	individuals.resize(npop, std::vector<T>(ndv));
	for (auto i = 0; i < npop; ++i)
	{
		individuals[i] = decision_variables;
	}
}