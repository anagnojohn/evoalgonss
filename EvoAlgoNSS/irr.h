#pragma once

#include "dependencies.h"
#include "bond.h"

template<typename T>
T irr(const T& r, const T& nominal_value, const std::vector<T>& cash_flows, const size_t& int_frequency)
{
	T frequency = static_cast<T>(int_frequency);
	const size_t& num_time_periods = cash_flows.size();
	T sum = 0.0;
	for (auto i = 0; i < num_time_periods; ++i)
	{
		sum = sum + cash_flows[i] / std::pow((1 + r / frequency), static_cast<T>(i + 1));
	}
	return sum + nominal_value / std::pow((1 + r / frequency), static_cast<T>(num_time_periods));
}

template<typename T>
T fitness_irr(const std::vector<T>& solution, const Bond<T>& bond)
{
	T sum_of_squares = 0;
	sum_of_squares = sum_of_squares + std::pow(bond.price - irr(solution[0], bond.nominal_value, bond.cash_flows, bond.frequency), 2);
	return sum_of_squares;
}