#pragma once

#include "dependencies.h"

template<typename T>
std::vector<T> construct_donor(const std::vector<std::vector<T>>& individuals)
{
	const auto& npop = individuals.size();
	const auto& ndv = individuals[0].size();
	std::vector<T> donor(ndv);
	std::vector<size_t> r_i;
	std::vector<size_t> indices;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	const T F_param = 0.4;
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
		donor[j] = individuals[r_i[0]][j] + F_param * (individuals[r_i[1]][j] - individuals[r_i[2]][j]);
	}
	return donor;
}

template<typename T>
std::vector<T> construct_trial(const std::vector<T>& target, const std::vector<T>& donor)
{
	const T Cr = 0.5;
	auto ndv = donor.size();
	std::vector<T> trial(ndv);
	std::random_device generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	std::vector<size_t> indices;
	std::random_device random_device;
	std::mt19937 engine{ random_device() };
	for (auto i = 0; i < ndv; ++i)
	{
		indices.push_back(i);
	}
	std::uniform_int_distribution<size_t> size_distribution(0, indices.size() - 1);
	for (auto j = 0; j < ndv; ++j)
	{
		T epsilon = distribution(generator);
		size_t jrand = indices[size_distribution(engine)];
		if (epsilon <= Cr || j == jrand)
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
std::vector<T> differential_evo(std::vector<std::vector<T>> individuals, F f, const T& tol, const T& opt, const size_t& gmax, const std::vector<T>& stdev)
{
	const size_t& npop = individuals.size();
	const size_t& ndv = individuals[0].size();
	std::random_device generator;
	std::uniform_real_distribution<> distribution(0.0, 1.0);
	boost::math::beta_distribution<> dist(1, 6);
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
	for (auto g = 0; g < gmax; ++g)
	{
		auto min_cost = f(individuals[0]);
		if (tol > std::abs(f(individuals[0]) - opt))
		{
			std::cout << "Found solution at iteration: " << g << "." << '\n';
			break;
		}
		for (auto i = 0; i < npop; ++i)
		{
			std::vector<T> donor = construct_donor(individuals);
			std::vector<T> trial = construct_trial(individuals[i], donor);
			if (f(trial) <= f(individuals[i]))
			{
				individuals[i] = trial;
			}
		}
	}
	std::sort(individuals.begin(), individuals.end(), comparator);
	std::cout << "Sorted Costs:" << '\n';
	for (const auto& p : individuals)
	{
	//	std::cout << p << " " << f(p) << '\n';
	}
	return individuals[0];
}