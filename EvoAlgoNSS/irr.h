#pragma once

#include <vector>
#include <cmath>

//! Internal Rate of Return (IRR)
namespace irr
{
	//! Calculates discount factors
	template<typename T>
	T compute_discount_factor(const T& r, const T& frequency, const T& period)
	{
		return 1 / std::pow((1 + r / frequency), period);
	}

	//! Constraints function for IRR
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

	//! Returns the present value of an investment
	template<typename T>
	T irr(const T& r, const T& nominal_value, const std::vector<T>& cash_flows, const T& frequency)
	{
		const size_t& num_time_periods = cash_flows.size();
		T sum = 0.0;
		for (auto i = 0; i < num_time_periods; ++i)
		{
			sum = sum + cash_flows[i] * compute_discount_factor(r, frequency, static_cast<T>(i + 1));
		}
		return sum + nominal_value * compute_discount_factor(r, frequency, static_cast<T>(num_time_periods));
	}

	//! Penalty function for IRR
	template<typename T>
	T penalty_irr(const std::vector<T>& solution)
	{
		const T& r = solution[0];
		const T C = 1000;
		T sum = 0;
		if (r < 0 || r > 1)
		{
			sum = sum + C * std::pow(std::abs(r), 2);
		}
		return sum;
	}

	//! This is the fitness function for finding the internal rate of return of a bond, in this case it is equal to its yield to maturity
	template<typename T>
	T fitness_irr(const std::vector<T>& solution, const T& price, const T& nominal_value, const std::vector<T>& cash_flows, const T& frequency)
	{
		T sum_of_squares = 0;
		sum_of_squares = sum_of_squares + std::pow(price - irr(solution[0], nominal_value, cash_flows, frequency), 2);
		sum_of_squares = sum_of_squares + penalty_irr(solution);
		return sum_of_squares;
	}
}