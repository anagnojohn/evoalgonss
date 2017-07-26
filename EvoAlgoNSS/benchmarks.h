#pragma once

#include <vector>
#include <tuple>
#include "bond.h"

template<typename T>
void print_results(std::tuple<std::vector<T>, T, size_t, double> results)
{
	std::cout << "Optimum solution: " << std::get<0>(results) << " Fitness Value: " << std::get<1>(results) << "\n";
	std::cout << "Solved at iteration: " << std::get<2>(results) << "\n";
	std::cout << "Elapsed time in seconds: " << std::get<3>(results) << "\n";
}

template<typename T, typename Sirr>
EAparams<T> set_init_nelson_param(std::vector<Bond<T>> & bonds, Sirr & irr_solver, const EAparams<T> & irr_params)
{
	for (auto i = 0; i < bonds.size(); ++i)
	{
		std::cout << "Processing bond: " << i + 1 << "\n";
		auto res = find_yield(bonds[i], irr_solver, irr_params);
		print_results(res);
		bonds[i].yield = std::get<0>(res)[0];
		bonds[i].duration = macaulay_duration(bonds[i].yield, bonds[i].cash_flows, bonds[i].nominal_value, bonds[i].frequency);
		std::cout << "Macaulay Duration: " << bonds[i].duration << "\n";
	}
	double b0 = 4.336 / 100;
	double b1 = 4.991 / 100 - b0;
	double b2 = 0;
	double b3 = 0;
	double tau1 = 1;
	double tau2 = tau1;
	const std::vector<T> stdev{ 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<T> decision_variables{ b0, b1, b2, b3, tau1, tau2 };
	const EAparams<T> nelson_params{ decision_variables, stdev };
	return nelson_params;
}

template<typename T, typename Sirr, typename Spricing>
void benchmarkbondpricing(std::vector<Bond<T>>& bonds, Sirr& irr_solver, Spricing& pricing_solver, const EAparams<T>& irr_params)
{
	const auto& pricing_params = set_init_nelson_param(bonds, irr_solver, irr_params);
	std::cout << "Solving bond pricing..." << "\n";
	print_results(bond_pricing(bonds, pricing_solver, pricing_params));
}

template<typename T, typename Sirr, typename Scurve>
void benchmarkcurvefitting(std::vector<Bond<T>>& bonds, Sirr& irr_solver, Scurve& curve_solver, const EAparams<T>& irr_params)
{
	const auto& curve_params = set_init_nelson_param(bonds, irr_solver, irr_params);
	std::cout << "Yield curve fitting..." << "\n";
	print_results(yield_curve_fitting(bonds, curve_solver, curve_params));
}