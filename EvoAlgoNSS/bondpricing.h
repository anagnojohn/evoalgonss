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
		sum = sum + coupon_value * std::exp(-svensson(solution, i + 1) * (i + 1));
	}
	return sum + nominal_value * std::exp(-svensson(solution, num_of_terms) * num_of_terms);
}

template<typename T>
T fitness_bond_pricing(const std::vector<T>& solution, const std::vector<Bond<T>>& bonds)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < bonds.size(); ++i)
	{
		sum_of_squares = sum_of_squares +
			std::pow(bonds[i].price - estimate_bond_pricing(solution, bonds[i].coupon_value, bonds[i].nominal_value, bonds[i].cash_flows.size()), 2)
			/ bonds[i].duration;
	}
	return sum_of_squares;
}

template<typename T, typename S>
std::vector<T> bond_pricing(std::vector< Bond<T> > bonds, S& solver)
{
	auto f = [&](const std::vector<double>& solution) { return fitness_bond_pricing(solution, bonds); };
	return solve(f, 0.0, solver);
}
