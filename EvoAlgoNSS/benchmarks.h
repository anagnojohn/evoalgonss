#pragma once

#include <vector>
#include <tuple>
#include "bond.h"
#include "differentialevo.h"

template<typename T, typename S>
std::tuple<std::vector<T>, T, size_t, double> yield_curve_fitting(const std::vector< Bond<T> >& bonds, const S& curve_solver, const Population<T>& curve_pop)
{
	assert(curve_pop.ndv == 6);
	for (const auto& p : bonds)
	{
		assert(p.yield > 0 && p.yield < 1);
		assert(p.duration > 0);
	}
	auto f = [&](const auto& solution) { return fitness_svensson(solution, bonds); };
	//auto solver = create_solver<T, decltype(f), S>(curve_solver, curve_pop);
	Solver<T, decltype(f), S> solver(curve_solver, curve_pop);
	return solver.solve(f, 0.0);
}

template<typename T, typename S>
std::tuple<std::vector<T>, T, size_t, double> bond_pricing(const std::vector< Bond<T> >& bonds, const S& pricing_solver, const Population<T>& pricing_pop)
{
	assert(curve_pop.ndv == 6);
	for (const auto& p : bonds)
	{
		assert(p.yield > 0 && p.yield < 1);
		assert(p.duration > 0);
	}
	auto f = [&](const auto& solution) { return fitness_bond_pricing(solution, bonds); };
	Solver<T, decltype(f), S> solver(pricing_solver, pricing_pop);
	return solver.solve(f, 0.0);
}

template<typename T, typename S>
std::tuple<std::vector<T>, T, size_t, double> find_yield(const Bond<T> bond, const S& irr_solver, const Population<T>& irr_pop)
{
	assert(irr_pop.ndv == 1);
	auto f = [&](const auto& solution) { return fitness_irr(solution, bond);};
	Solver<T, decltype(f), S> solver(irr_solver, irr_pop);
	return solver.solve(f, 0.0);
}

template<typename T>
void print_results(std::tuple<std::vector<T>, T, size_t, double> results)
{
	std::cout << "Optimum solution: " << std::get<0>(results) << " Fitness Value: " << std::get<1>(results) << "\n";
	std::cout << "Solved at iteration: " << std::get<2>(results) << "\n";
	std::cout << "Elapsed time in seconds: " << std::get<3>(results) << "\n";
}

template<typename T, typename Sirr>
Population<T> set_init_nelson_param(std::vector<Bond<T>> & bonds, const Sirr& irr_solver, const Population<T>& irr_pop, const size_t& npop)
{
	for (auto i = 0; i < bonds.size(); ++i)
	{
		std::cout << "Processing bond: " << i + 1 << "\n";
		auto res = find_yield(bonds[i], irr_solver, irr_pop);
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
	const Population<T> nelson_params{ decision_variables, stdev, npop};
	return nelson_params;
}

template<typename T, typename Sirr, typename Spricing>
void benchmarkbondpricing(std::vector<Bond<T>>& bonds, const Sirr& irr_solver, const Spricing& pricing_solver)
{
	Population<T> irr_pop{ { 0.2 },{ 0.7 }, 200 };
	const auto& pricing_pop = set_init_nelson_param(bonds, irr_solver, irr_pop, 200);
	std::cout << "Solving bond pricing..." << "\n";
	print_results(bond_pricing(bonds, pricing_solver, pricing_pop));
}

template<typename T, typename Sirr, typename Scurve>
void benchmarkcurvefitting(std::vector<Bond<T>>& bonds, const Sirr& irr_solver, const Scurve& curve_solver)
{
	Population<T> irr_pop{ { 0.2 },{ 0.7 }, 200 };
	const auto& curve_pop= set_init_nelson_param(bonds, irr_solver, irr_pop, 200);
	std::cout << "Yield curve fitting..." << "\n";
	print_results(yield_curve_fitting(bonds, curve_solver, curve_pop));
}