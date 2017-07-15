#pragma once

#include "dependencies.h"

template<typename T, typename F>
class DifferentialEvo
{
public:
	DifferentialEvo(std::vector<std::vector<T>> i_individuals, const T& i_tol, const T& i_opt, const size_t& i_gmax, const T& i_cr, const T& i_f_param, const std::vector<T>& i_stdev)
		: individuals(i_individuals), cr(i_cr), f_param(i_f_param), tol(i_tol), opt(i_opt), gmax(i_gmax), stdev(i_stdev), npop(individuals.size()), ndv(individuals[0].size())
	{
		std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
		distribution = i_distribution;
	}
	std::vector<T> differential_evo(F f);
private:
	std::vector<std::vector<T>> individuals;
	const T cr; // 0.5;
	const T f_param; // 0.4;
	const T tol;
	const T opt;
	const std::vector<T> stdev;
	const size_t gmax;
	size_t npop;
	const size_t ndv;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	std::random_device generator;
	std::uniform_real_distribution<T> distribution;
	std::vector<T> construct_donor();
	std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor);
};

template<typename T, typename F>
std::vector<T> DifferentialEvo<T, F>::construct_donor()
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

template<typename T, typename F>
std::vector<T> DifferentialEvo<T, F>::construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
{
	
	auto ndv = donor.size();
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
std::vector<T> DifferentialEvo<T, F>::differential_evo(F f)
{
	init_epsilon(individuals, stdev);
	//std::cout << "Individuals:" << '\n';
	//for (const auto& p : individuals)
	//{
	//	std::cout << p << " " << f(p) << '\n';
	//}
	auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
	{
		return f(l) < f(r);
	};
	T min_cost = f(individuals[0]);
	for (const auto& p : individuals)
	{
		if (min_cost > f(p))
		{
			min_cost = f(p);
		}
	}
	for (auto g = 0; g < gmax; ++g)
	{
		for (const auto& p : individuals)
		{
			if (min_cost > f(p))
			{
				min_cost = f(p);
			}
		}
		//std::cout << min_cost << "\n";
		if (tol > std::abs(min_cost - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		for (auto i = 0; i < npop; ++i)
		{
			std::vector<T> donor = construct_donor();
			std::vector<T> trial = construct_trial(individuals[i], donor);
			if (f(trial) <= f(individuals[i]))
			{
				individuals[i] = trial;
			}
		}
	}
	T error = std::abs(min_cost - opt);
	std::sort(individuals.begin(), individuals.end(), comparator);
	//std::cout << "Sorted Costs:" << '\n';
	//for (const auto& p : individuals)
	//{
	//	std::cout << p << " " << f(p) << '\n';
	//}
	return individuals[0];
}