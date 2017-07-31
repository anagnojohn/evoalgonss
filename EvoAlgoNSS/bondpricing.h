#pragma once

#include "dependencies.h"
#include "svensson.h"
#include "irr.h"
#include "bond.h"

template<typename T>
T estimate_bond_pricing(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, size_t num_of_terms)
{
	T sum = 0;
	for (auto i = 0; i < num_of_terms; ++i)
	{
		sum = sum + coupon_value * std::exp(-svensson(solution, static_cast<T>(i + 1)) * static_cast<T>(i + 1));
	}
	return sum + nominal_value * std::exp(-svensson(solution, static_cast<T>(num_of_terms)) * static_cast<T>(num_of_terms));
}

template<typename T>
T fitness_bond_pricing(const std::vector<T>& solution, const std::vector<Bond<T>>& bonds)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < bonds.size(); ++i)
	{
		sum_of_squares = sum_of_squares +
			std::pow(bonds[i].price - estimate_bond_pricing(solution, bonds[i].coupon_value, bonds[i].nominal_value, bonds[i].cash_flows.size()), 2)
			/ std::sqrt(bonds[i].duration);
	}
	return sum_of_squares;
}