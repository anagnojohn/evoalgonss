#pragma once

#include <vector>
#include <assert.h>

template<typename T, typename F>
class EA_base
{
public:
	EA_base(const size_t& i_npop, const T& i_tol, const size_t& i_iter_max) : tol{i_tol}, iter_max(i_iter_max), npop{i_npop}
	{
		assert(npop > 0);
		assert(tol > 0);
		assert(iter_max > 0);
		//if (individuals.size() > 0)
		//{
		//	individuals.erase(individuals.begin(), individuals.end());
		//}
		//init_individuals();
	}
	// Getters for parameters
	size_t get_npop() const { return npop; ; };
	std::vector<std::vector<T>> get_individuals() const { return individuals; };
	T get_tol() const { return tol; };
	size_t get_iter_max() const { return iter_max; };
	// Setters for parameters
	void set_npop(const size_t& npop) { assert(npop > 0); this->npop = npop; create_individuals(); };
	void set_tol(const T& tol) { assert(tol > 0); this->tol = tol; };
	void set_iter_max(const size_t& iter_max) { assert(iter_max > 0); this->iter_max = iter_max; };
	std::vector<std::vector<T>> init_individuals(const std::vector<T>& decision_variables, const std::vector<T>& stdev);
protected:
	// Tolerance
	T tol;
	// Number of maximum iterations
	size_t iter_max;
	// Size of the population
	size_t npop;
	// Population
	std::vector<std::vector<T>> individuals;
};

template<typename T, typename F>
std::vector<std::vector<T>> EA_base<T, F>::init_individuals(const std::vector<T>& decision_variables, const std::vector<T>& stdev)
{
	std::random_device generator;
	T epsilon;
	const auto& ndv = decision_variables.size();
	std::vector<std::vector<T>> individuals;
	individuals.resize(npop, std::vector<T>(ndv));
	for (auto& p : individuals)
	{
		p = decision_variables;
	}
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			std::normal_distribution<T> ndistribution(0, stdev[j]);
			epsilon = ndistribution(generator);
			individuals[i][j] = individuals[i][j] + epsilon;
		}
	}
	return individuals;
}

template<typename T>
class EAparams
{
public:
	EAparams(const std::vector<T>& i_decision_variables, const std::vector<T>& i_stdev)
	: decision_variables{ i_decision_variables }, stdev{i_stdev}
	{
		assert(decision_variables.size() > 0);
		assert(decision_variables.size() == stdev.size());
		for (const auto& p : stdev)
		{
			assert(p > 0);
		}
		this->decision_variables = decision_variables;
		this->ndv = decision_variables.size();
		this->stdev = stdev;
	}
	// Getters for parameters
	size_t get_ndv() const { return ndv; };
	std::vector<T> get_stdev() const { return stdev; };
	std::vector<T> get_decision_variables() const { return decision_variables; };
	// Setters for parameters
	void set_stdev(const std::vector<T>& stdev)
	{
		assert(decision_variables.size() == stdev.size());
		for (const auto& p : stdev)
		{
			assert(p > 0);
		}
		this->stdev = stdev;
	};
	void set_decision_variables(const std::vector<T>& decision_variables)
	{
		this->decision_variables = decision_variables; ndv = decision_variables.size();
	};
private:
	// Standard deviation of the decision variables
	std::vector<T> stdev;
	// Decision Variables
	std::vector<T> decision_variables;
	// Number of decision variables
	size_t ndv;
};