#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

template<typename T>
struct Interest_Rate
{
	Interest_Rate(const T& i_period, const T& i_rate) : period{ i_period }, rate{ i_rate } {};
	const T period;
	const T rate;
};

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

template<typename T>
T fitness_yield_curve_fitting(const std::vector<T>& solution, const std::vector<Interest_Rate<T>>& ir_vec)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < ir_vec.size(); ++i)
	{
		T estimate = svensson(solution, ir_vec[i].period);
		sum_of_squares = sum_of_squares + std::pow(ir_vec[i].rate - estimate, 2);
	}
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
	const T C = 1000;
	if (b0 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0), 2);
	}
	if (b0 + b1 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0 + b1), 2);
	}
	if (tau1 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	if (tau2 < 0)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	return sum_of_squares;
}