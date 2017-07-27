#pragma once

#include "dependencies.h"
#include "bond.h"

template<typename T>
T svensson(const std::vector<T>& solution, const T& m)
{
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
		return 100;
	}
}

template<typename T>
T fitness_svensson(const std::vector<T>& solution, const std::vector<Bond<T>>& bonds)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < bonds.size(); ++i)
	{
		sum_of_squares = sum_of_squares + std::pow(bonds[i].yield - svensson(solution, bonds[i].duration), 2) / std::sqrt(bonds[i].duration);
	}
	return sum_of_squares;
}