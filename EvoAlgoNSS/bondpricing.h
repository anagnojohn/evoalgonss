#pragma once

#include "dependencies.h"
#include "svensson.h"
#include "irr.h"

template<typename T>
T macaulay_duration(T bond_price)
{
	return;
}

template<typename T>
T estimate_bond_price(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, size_t num_of_terms, T term, T maturity)
{
	T sum = 0;
	for (auto i = 0; i < terms; ++i)
	{
		sum = sum + coupon_value * *std::exp(-svensson(solution, term) * t);
	}
	return sum + nominal_value * std::exp(-svensson(solution, maturity * m);
}

template<typename T>
T fitness_bond_price(const std::vector<T>& solution, const std::vector<T>& bond_prices, F svensson, const T& coupon_value, const T& nominal_value, size_t num_of_terms, T term, const std::vector<T>& maturity)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < bond_prices.size(); ++i)
	{
		sum_of_squares = sum_of_squares +
			std::pow(bond_prices[i] - estimate_bond_prices(solution, coupon_value, nominal_value, num_of_terms, term, maturity[i])
				* 1 / (1 / macaulay_duration(bond_prices[i]))
	}
	return sum_of_squares;
}