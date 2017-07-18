#pragma once

#include "dependencies.h"

template<typename T>
class DifferentialEvo
{
public:
	DifferentialEvo(const std::vector<T>& i_decision_variables, const size_t& i_npop, const T& i_tol, const size_t& i_gmax, const T& i_cr, const T& i_f_param, const std::vector<T>& i_stdev)
		: cr(i_cr), f_param(i_f_param), tol(i_tol), gmax(i_gmax), stdev(i_stdev), npop(i_npop), ndv(i_decision_variables.size())
	{
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
		individuals = create_individuals(npop, i_decision_variables);
		init_epsilon(individuals, stdev);
	}
	std::vector<std::vector<T>> individuals;
	const T tol;
	std::vector<T> construct_donor();
	std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
	const size_t gmax;
private:
	const T cr; // 0.5;
	const T f_param; // 0.4;
	const std::vector<T> stdev;
	size_t npop;
	const size_t ndv;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
};

template<typename T>
std::vector<T> DifferentialEvo<T>::construct_donor()
{
	std::vector<T> donor(ndv);
	std::vector<size_t> r_i;
	std::vector<size_t> indices;
	for (auto i = 0; i < npop; ++i)
	{
		indices.push_back(i);
	}
	std::uniform_int_distribution<size_t> distribution(0, indices.size() - 1);
	while (r_i.size() < 4)
	{
		r_i.push_back(indices[distribution(engine)]);
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
std::vector<T> DifferentialEvo<T>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
{
	std::vector<T> trial(ndv);
	std::vector<size_t> indices;
	for (auto i = 0; i < ndv; ++i)
	{
		indices.push_back(i);
	}
	std::uniform_int_distribution<size_t> size_distribution(0, indices.size() - 1);
	for (auto j = 0; j < ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = indices[size_distribution(engine)];
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

template<typename T, typename F>
std::vector<T> solve(F f, const T& opt, DifferentialEvo<T>& de)
{
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	std::vector<T> min_cost = de.individuals[0];
	for (const auto& p : de.individuals)
	{
		if (f(min_cost) > f(p))
		{
			min_cost = p;
		}
	}
	for (auto g = 0; g < de.gmax; ++g)
	{
		if (de.tol > std::abs(f(min_cost) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		for (auto i = 0; i < de.individuals.size(); ++i)
		{
			std::vector<T> donor = de.construct_donor();
			std::vector<T> trial = de.construct_trial(de.individuals[i], donor);
			if (f(trial) <= f(de.individuals[i]))
			{
				de.individuals[i] = trial;
			}
		}
		for (const auto& p : de.individuals)
		{
			if (f(min_cost) > f(p))
			{
				min_cost = p;
			}
		}
	}
	T error = std::abs(f(min_cost) - opt);
	return min_cost;
}