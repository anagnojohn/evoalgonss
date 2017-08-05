#pragma once

#include "dependencies.h"
#include "bond.h"

template<typename T>
bool constraints_svensson(const std::vector<T>& solution)
{
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
	if (b0 > 0 && b0 + b1 > 0 && tau1 > 0 && tau2 > 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}
// Spot interest rate at term m
template<typename T>
T svensson(const std::vector<T>& solution, const T& m)
{
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
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

template<typename T>
void penalty_svensson(const std::vector<T>& solution, T& sum_of_squares)
{
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
	const T C = 1000;
	if (b0 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0), 2);
	}
	if (b0 + b1 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0 + b1), 2);
	}
	if (tau1 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	if (tau2 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	return sum_of_squares;
}