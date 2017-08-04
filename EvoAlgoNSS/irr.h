#pragma once

#include "dependencies.h"
#include "bond.h"

template<typename T>
bool constraints_irr(const std::vector<T>& solution)
{
	const T& r = solution[0];
	if (r > 0 && r < 1)
	{
		return true;
	}
	else
	{
		return false;
	}
}


template<typename T>
T irr(const T& r, const T& nominal_value, const std::vector<T>& cash_flows, const T& frequency)
{
	const size_t& num_time_periods = cash_flows.size();
	T sum = 0.0;
	for (auto i = 0; i < num_time_periods; ++i)
	{
		sum = sum + cash_flows[i] / std::pow((1 + r / frequency), static_cast<T>(i + 1));
	}
	return sum + nominal_value / std::pow((1 + r / frequency), static_cast<T>(num_time_periods));
}

template<typename T>
T find_bond_price(const T& ytm, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods)
{
	const size_t& num_time_periods = time_periods.size();
	T sum = 0.0;
	T pv_cash_flow = 0.0;
	T discount_factor = 0.0;
	for (auto i = 0; i < num_time_periods; ++i)
	{
		discount_factor = std::exp(-ytm * time_periods[i]);
		pv_cash_flow = pv_cash_flow + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + nominal_value * discount_factor;
	return pv_cash_flow;
}

template<typename T>
T fitness_irr(const std::vector<T>& solution, const Bond<T>& bond)
{
	T sum_of_squares = 0;
	sum_of_squares = sum_of_squares + std::pow(bond.price - irr(solution[0], bond.nominal_value, bond.cash_flows, bond.frequency), 2);
	return sum_of_squares;
}