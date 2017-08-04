#pragma once

#include <vector>
#include <tuple>
#include "bond.h"
#include "svensson.h"
#include "bondpricing.h"
#include "irr.h"
#include "yield_curve_fitting.h"
#include "geneticalgo.h"
#include "local_best_pso_variants.h"
#include "differentialevo.h"

template<typename T>
std::vector<T> set_test_bond_data(std::vector<Bond<T>>& bonds)
{
	bonds[0].yield = 0.054308895;
	bonds[1].yield = 0.090624152;
	bonds[2].yield = 0.030896968;
	bonds[3].yield = 0.006625537;
	bonds[4].yield = 0.07972484;
	bonds[5].yield = 0.03366204;
	bonds[6].yield = 0.039963969;
	bonds[7].yield = 0.070339142;
	bonds[0].duration = 2.944988754;
	bonds[1].duration = 4.966178711;
	bonds[2].duration = 6.279674883;
	bonds[3].duration = 9.474865358;
	bonds[4].duration = 8.416273259;
	bonds[5].duration = 14.93089635;
	bonds[6].duration = 15.38446779;
	bonds[7].duration = 12.22184684;
	double b0 = bonds[0].yield;
	double b1 = bonds[7].yield - b0;
	double b2 = 0;
	double b3 = 0;
	double tau1 = bonds[5].duration;
	double tau2 = tau1;
	const std::vector<T> decision_variables{ b0, b1, b2, b3, tau1, tau2 };
	return decision_variables;
}

template<typename T>
void print_results(std::tuple<std::vector<T>, T, size_t, double> results)
{
	std::cout << "Optimum solution: " << std::get<0>(results) << " Fitness Value: " << std::get<1>(results) << "\n";
	std::cout << "Solved at iteration: " << std::get<2>(results) << "\n";
	std::cout << "Elapsed time in seconds: " << std::get<3>(results) << "\n";
}

template<typename T, typename S>
std::vector<T> set_init_nelson_param(std::vector<Bond<T>> & bonds, const S& irr_solver)
{
	assert(irr_solver.ndv == 1);
	for (auto i = 0; i < bonds.size(); ++i)
	{
		std::cout << "Processing bond: " << i + 1 << "\n";
		auto f = [&](const auto& solution) { return fitness_irr(solution, bonds[i]);};
		Solver<T, decltype(f), S> solver{ irr_solver };
		auto res = solver.solve(f, 0.0);
		print_results(res);
		bonds[i].yield = std::get<0>(res)[0];
		bonds[i].duration = macaulay_duration(bonds[i].yield, bonds[i].cash_flows, bonds[i].nominal_value, bonds[i].frequency);
		std::cout << "Macaulay Duration: " << bonds[i].duration << "\n";
	}
	double b0 = bonds[0].yield;
	double b1 = bonds[7].yield - b0;
	double b2 = 0;
	double b3 = 0;
	double tau1 = bonds[5].duration;
	double tau2 = tau1;
	const std::vector<T> decision_variables{ b0, b1, b2, b3, tau1, tau2 };
	return decision_variables;
}

template<typename T, typename S>
void benchmarkbondpricing(std::vector<Bond<T>>& bonds,  const S& pricing_solver)
{
	assert(pricing_solver.ndv == 6);
	for (const auto& p : bonds)
	{
		assert(p.yield > 0 && p.yield < 1);
		assert(p.duration > 0);
	}
	auto f = [&](const auto& solution) { return fitness_bond_pricing(solution, bonds); };
	Solver<T, decltype(f), S> solver{ pricing_solver };
	std::cout << "Solving bond pricing using bond prices..." << "\n";
	auto res = solver.solve(f, 0.0);
	print_results(res);
	for (const auto& p : bonds)
	{
		std::cout << "Estimated yield: " << svensson(std::get<0>(res), p.duration) << " Actual Yield: " << p.yield << "\n";
	}
	for (const auto& p : bonds)
	{
		std::cout << "Estimated price: " << estimate_bond_pricing(std::get<0>(res), p.coupon_value, p.nominal_value, p.time_periods) << " Actual Price: " << p.price << "\n";
	};
}

template<typename T, typename S>
void benchmarkbondpricing_yields(std::vector<Bond<T>>& bonds, const S& pricing_solver)
{
	assert(pricing_solver.ndv == 6);
	for (const auto& p : bonds)
	{
		assert(p.yield > 0 && p.yield < 1);
		assert(p.duration > 0);
	}
	auto f = [&](const auto& solution) { return fitness_bond_pricing_yields(solution, bonds); };
	Solver<T, decltype(f), S> solver{ pricing_solver };
	std::cout << "Solving bond pricing using bond yields..." << "\n";
	auto res = solver.solve(f, 0.0);
	print_results(res);
	for (const auto& p : bonds)
	{
		std::cout << "Estimated yield: " << svensson(std::get<0>(res), p.duration) << " Actual Yield: " << p.yield << "\n";
	}
	for (const auto& p : bonds)
	{
		std::cout << "Estimated price: " << estimate_bond_pricing(std::get<0>(res), p.coupon_value, p.nominal_value, p.time_periods) << " Actual Price: " << p.price << "\n";
	}
}

template<typename T, typename S>
void benchmarkyieldcurvefitting(const std::vector<Interest_Rate<T>>& ir_vec, const S& curve_solver)
{
	assert(curve_solver.ndv == 6);
	auto f = [&](const auto& solution) { return fitness_yield_curve_fitting(solution, ir_vec); };
	Solver<T, decltype(f), S> solver{ curve_solver };
	std::cout << "Yield Curve fitting." << "\n";
	auto res = solver.solve(f, 0.0);
	print_results(res);
	for (const auto& p : ir_vec)
	{
		std::cout << "Estimated interest rates: " << svensson(std::get<0>(res), p.period) << " Actual interest rates: " << p.rate << "\n";
	}
}