#pragma once

#include <vector>
#include <cmath>

//! Internal Rate of Return (IRR)
namespace irr
{
	//! Calculates discount factors
	template<typename T>
	T compute_discount_factor(const T& r, const T& period, const DF_type& df_type)
	{
		switch (df_type)
		{
		case (DF_type::frac): return 1 / std::pow((1 + r), period);
		case (DF_type::exp): return std::exp(-r * period);
		}
	}

	//! Constraints function for IRR
	template<typename T>
	bool constraints_irr(const std::vector<T>& solution, const Constraints_type& constraints_type)
	{
		switch (constraints_type)
		{
		case(Constraints_type::normal):
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
		case(Constraints_type::tight): return true;
		case(Constraints_type::none): return true;
		}
	}

	//! Returns the present value of an investment
	template<typename T>
	T compute_pv(const T& r, const T& nominal_value, const std::vector<T>& cash_flows, const std::vector<T>& time_periods, const DF_type& df_type)
	{
		assert(cash_flows.size() == time_periods.size());
		const size_t& num_time_periods =time_periods.size();
		T sum = 0.0;
		for (size_t i = 0; i < num_time_periods; ++i)
		{
			sum = sum + cash_flows[i] * compute_discount_factor(r, time_periods[i], df_type);
		}
		return sum + nominal_value * compute_discount_factor(r, time_periods.back(), df_type);
	}

	//! Penalty function for IRR
	template<typename T>
	T penalty_irr(const T& r)
	{
		T sum = 0;
		const T C = 1000;
		if (r < 0 || r > 1)
		{
			sum = sum + C * std::pow(std::abs(r), 2);
		}
		return sum;
	}

	//! This is the fitness function for finding the internal rate of return of a bond, in this case it is equal to its yield to maturity
	template<typename T>
	T fitness_irr(const std::vector<T>& solution, const T& price, const T& nominal_value, const std::vector<T>& cash_flows, const std::vector<T>& time_periods, 
		const DF_type& df_type, const bool& use_penalty_method)
	{
		T sum_of_squares = 0;
		sum_of_squares = sum_of_squares + std::pow(price - compute_pv(solution[0], nominal_value, cash_flows, time_periods, df_type), 2);
		if (use_penalty_method)
		{
			return sum_of_squares + penalty_irr(solution[0]);
		}
		else
		{
			return sum_of_squares;
		}
	}
}