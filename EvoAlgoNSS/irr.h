#pragma once

#include <vector>
#include <cmath>
#include "utilities.h"

using namespace utilities;

//! Internal Rate of Return (IRR) namespace
namespace irr
{
	/** \fn compute_discount_factor(const T& r, const T& period, const DF_type& df_type)
	*  \brief Calculates discount factors
	*  \param r Rate
	*  \param period The period the rate was recorded
	*  \param df_type The method used to calculate the discount factor
	*  \return The discount factor
	*/
	template<typename T>
	T compute_discount_factor(const T& r, const T& period, const DF_type& df_type)
	{
		switch (df_type)
		{
		case (DF_type::frac): return 1 / std::pow((1 + r), period);
		case (DF_type::exp): return std::exp(-r * period);
            default: return std::exp(-r * period);
		}
	}

	/** \fn constraints_irr(const std::vector<T>& solution, const Constraints_type& constraints_type)
	*  \brief Constraints function for Internal Rate of Return
	*  \param solution Internal Rate of Return candindate solution
	*  \param constraints_type Type of constraints used
	*  \return True if constraints are satisfied, false otherwise
	*/
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
            default: return true;
		}
	}

	/** \fn compute_pv(const T& r, const T& nominal_value, const std::vector<T>& cash_flows, const std::vector<T>& time_periods, const DF_type& df_type)
	*  \brief Returns the present value of an investment
	*  \param r Internal Rate of Return
	*  \param nominal_value The nominal value of the investment
	*  \param cash_flows The cash flows of the investment
	*  \param time_periods The time periods that correspond to the cash flows of the investment
	*  \param df_type The method used to calculate the discount factor
	*  \return The present value of the investment
	*/
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

	/** \fn penalty_irr(const T& r)
	*  \brief Penalty function for IRR
	*  \param r Candidate solution for the Internal Rate of Return
	*  \return A penalty value, if constraints are not satisfied
	*/
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

	/** \fn fitness_irr(const std::vector<T>& solution, const T& price, const T& nominal_value, const std::vector<T>& cash_flows, const std::vector<T>& time_periods, 
		const DF_type& df_type, const bool& use_penalty_method)
	*  \brief This is the fitness function for finding the internal rate of return of a bond, in this case it is equal to its yield to maturity
	*  \param solution Internal Rate of Return candindate solution
	*  \param price The present value of the investment
	*  \param nominal_value The nominal value of the investment
	*  \param cash_flows The cash flows of the investment
	*  \param time_periods The time periods that correspond to the cash flows of the investment
	*  \param df_type The method used to calculate the discount factor
	*  \param use_penalty_method Whether to use the penalty method defined for IRR or not
	*  \return The fitness cost of IRR
	*/
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