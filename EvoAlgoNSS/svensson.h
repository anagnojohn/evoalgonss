#pragma once

#include "dependencies.h"
#include "bond.h"

template<typename T>
T svensson(const std::vector<T>& solution, const size_t& int_m)
{
	T m = static_cast<T>(int_m);
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
	if (b0 > 0 && b0 + b1 > 0 && tau1 > 0 && tau2 > 0)
	{
		if (m == 0)
		{
			return b0 + b1;
		}
		else
		{
			T result = b0 + (b1 + b2) * (tau1 / m) * (1 - std::exp(-m / tau1))
				- b2 * std::exp(-m / tau1) + b3 * (tau2 / m) * (1 - std::exp(-m / tau2))
				- b3 * std::exp(-m / tau2);
			return result;
		}
	}
	else
	{
		return 1000;
	}
}

template<typename T>
T fitness_svensson(const std::vector<T>& solution, const std::vector<Bond<T>>& bonds)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < bonds.size(); ++i)
	{
		sum_of_squares = sum_of_squares + std::pow(bonds[i].yield - svensson(solution, bonds[i].cash_flows.size()), 2);
	}
	return sum_of_squares;
}

template<typename T, typename S>
std::vector<T> yield_curve_fitting(std::vector< Bond<T> > bonds, S& solver)
{
	auto f = [&](const std::vector<double>& solution) { return fitness_svensson<double>(solution, bonds); };
	return solve(f, 0.0, solver);
}