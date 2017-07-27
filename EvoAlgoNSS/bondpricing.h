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

/*
template<typename T>
std::tuple<std::vector<T>, T, size_t, double> bond_pricing(std::vector< Bond<T> > bonds)
{
	for (const auto& p : bonds)
	{
		assert(p.yield > 0 && p.yield < 1);
		assert(p.duration > 0);
	}
	auto f = [&](const auto& solution) { return fitness_bond_pricing(solution, bonds); };
	S solver{ ...,f };
	solver.solve();
	return solve(f, 0.0, solver, ea);
}
*/
