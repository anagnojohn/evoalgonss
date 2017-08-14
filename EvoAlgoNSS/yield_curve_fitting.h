#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "svensson.h"

using namespace nss;

//! Struct for interest rates
template<typename T>
struct Interest_Rate
{
	//! Constructor
	Interest_Rate(const T& i_period, const T& i_rate) : period{ i_period }, rate{ i_rate } {};
	//! Period is the time that matches the rate
	const T period;
	//! Interest rate
	const T rate;
};

//! Reads the interest rates and periods from file and constructs a vector of interest rate structs
template<typename T>
std::vector<Interest_Rate<T>> read_ir_from_file(const std::string & filename)
{
	std::vector<Interest_Rate<T>> ir_vec;
	std::ifstream input(filename);
	for (std::string line; getline(input, line); )
	{
		T period;
		T rate;
		std::istringstream stream(line);
		stream >> period >> rate;
		const Interest_Rate<T> ir{ period, rate / 100.0 };
		ir_vec.push_back(ir);
	}
	return ir_vec;
}

//! A class for the bond pricing problem
template<typename T>
class InterestRate_Helper
{
public:
	//! Constructor
	InterestRate_Helper(const std::vector<Interest_Rate<T>>& i_ir_vec) : ir_vec{ i_ir_vec } {};
	//! Yield Curve Fitting using interest rates and recorded periods
	template<typename S>
	void yieldcurve_fitting(const S& solver)
	{
		assert(solver.ndv == 6);
		auto f = [&](const auto& solution) { return fitness_yield_curve_fitting(solution); };
		auto c = [&](const auto& solution) { return constraints_svensson(solution); };
		std::cout << "Yield Curve fitting." << "\n";
		auto res = solve(f, c, solver);
		for (const auto& p : ir_vec)
		{
			std::cout << "Estimated interest rates: " << svensson(res, p.period) << " Actual interest rates: " << p.rate << "\n";
		}
	};
private:
	//! Vector of interest rates
	std::vector<Interest_Rate<T>> ir_vec;
	//! This is the fitness function for yield-curve fitting using Interest Rates
	T fitness_yield_curve_fitting(const std::vector<T>& solution)
	{
		//! The sum of squares of errors betwwen the actual rates and the rates computed by svensson are used
		T sum_of_squares = 0;
		for (size_t i = 0; i < ir_vec.size(); ++i)
		{
			T estimate = svensson(solution, ir_vec[i].period);
			sum_of_squares = sum_of_squares + std::pow(ir_vec[i].rate - estimate, 2);
		}
		sum_of_squares = sum_of_squares + penalty_svensson(solution);
		return sum_of_squares;
	};
};