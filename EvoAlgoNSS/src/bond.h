/** \file bond.h
* \author Ioannis Anagnostopoulos
* \brief Classes and functions for bonds and their internal rate of return
*/

#pragma once

#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>
#include <locale>
#include <assert.h>
#include "date.h"
#include "irr.h"
#include "ealgorithm_base.h"



//! Bond Class and Utilities
namespace bond
{
	using namespace ea;
	using namespace irr;
	template<typename T>
	class BondHelper;

	/*! \class Bond
	*  \brief Bond Class definition
	*/
	template<typename T>
	class Bond
	{
		//! Friend class BondHelper
		friend class BondHelper<T>;
	public:
		/** \fn Bond(const T& i_coupon_percentage, const T& i_price, const T& i_nominal_value, const T& i_frequency,
			std::string& i_settlement_date, const std::string& i_maturity_date)
		*  \brief Constructor
		*  \param i_coupon_percentage The coupon rate in %
		*  \param i_price The price of the bond
		*  \param i_nominal_value The nominal value of the bond
		*  \param i_frequency The frequency of coupon payments per year
		*  \param i_settlement_date The date the bond was bought
		*  \param i_maturity_date The date the bond expires
		*  \return A Bond<T> object
		*/
		Bond(const T& i_coupon_percentage, const T& i_price, const T& i_nominal_value, const T& i_frequency,
			std::string& i_settlement_date, const std::string& i_maturity_date) :
			coupon_percentage{ i_coupon_percentage },
			price{ i_price },
			nominal_value{ i_nominal_value },
			frequency{ i_frequency },
			coupon_value{ coupon_percentage * nominal_value / frequency },
			yield{ 0 },
			duration{ 0 }
		{
			assert(price > 0);
			assert(coupon_percentage > 0 && coupon_percentage < 1);
			assert(nominal_value > 0);
			assert(frequency > 0);
			std::tm t1 {};
			std::tm t2 {};
            std::stringstream s1;
            std::stringstream s2;
            s1 << i_settlement_date;
            s2 << i_maturity_date;
			s1 >> std::get_time(&t1, "%Y-%m-%d");
			s2 >> std::get_time(&t2, "%Y-%m-%d");
			settlement_date = date::year(t1.tm_year) / (t1.tm_mon+1) / t1.tm_mday;
			maturity_date = date::year(t2.tm_year) / (t2.tm_mon+1) / t2.tm_mday;
			cash_flows = compute_cash_flows();
			time_periods.resize(cash_flows.size());
			for (size_t i = 0; i < time_periods.size(); ++i)
			{
				time_periods[i] = static_cast<T>(i + 1) / frequency;
			}
		}
		/** \fn compute_yield(const T& i_price, const S& solver, const DF_type& df_type) const
		*  \brief Calculates the yield-to-maturity using the supplied solver
		*  \param i_price The price of the bond
		*  \param solver The parameter structure of the solver that is going to be used to estimate the yield of maturity
		*  \param df_type The type of discount factor method
		*  \return The yield-to-maturity of the bond
		*/
		template<typename S> T compute_yield(const T& i_price, const S& solver, const DF_type& df_type) const;
		/** \fn compute_yield(const T& i_price, const S& solver, const DF_type& df_type, const std::string& bonds_identifier) const
		*  \brief Calculates the yield-to-maturity using the supplied solver and passes the bond identifier to the solver
		*  \param i_price The price of the bond
		*  \param solver The parameter structure of the solver that is going to be used to estimate the yield of maturity
		*  \param df_type The type of discount factor method
		*  \param bonds_identifier An identifier for the bond in std::string form
		*  \return The yield-to-maturity of the bond
		*/
		template<typename S> T compute_yield(const T& i_price, const S& solver, const DF_type& df_type, const std::string& bonds_identifier) const;
		/** \fn compute_macaulay_duration(const DF_type& df_type)
		*  \brief Calculates the Macaulay duration of the bond
		*  \param df_type The type of discount factor method
		*  \return The Macaulay Duration of the bond
		*/
		T compute_macaulay_duration(const DF_type& df_type) const;
	private:
		/** \brief Bond's annual coupon rate */
		const T coupon_percentage;
		/** \brief Bond's price */
		const T price;
		/** \brief Bond's face value */
		const T nominal_value;
		/** \brief Bond's coupon payment frequency */
		const T frequency;
		/** \brief This is the annual coupon divided by the frequency */
		const T coupon_value;
		/** \brief Coupon payment periods */
		std::vector<T> time_periods;
		/** \brief A vector with all the coupon payments corresponding to time periods */
		std::vector<T> cash_flows;
		/** \brief Settlement date of the bond */
		date::sys_days settlement_date;
		/** \brief Maturity date of the bond */
		date::sys_days maturity_date;
		/** \brief Yield-to-maturity of the bond */
		T yield;
		/** \brief Macaulay duration of the bond */
		T duration;
		/** \fn compute_cash_flows()
		*  \brief Calculate the cash flows of the bond
		*  \return The cash flows of the bonds (coupon payments)
		*/
		std::vector<T> compute_cash_flows();
	};

	template<typename T>
	std::vector<T> Bond<T>::compute_cash_flows()
	{
		assert(settlement_date < maturity_date);
		const T number_of_days_coupon = 365 / frequency;
		const auto days_difference = (maturity_date - settlement_date).count();
		const auto time_periods = static_cast<T>(days_difference) / number_of_days_coupon;
		const size_t num_time_periods = static_cast<size_t>(std::ceil(time_periods));
		std::vector<T> cash_flows(num_time_periods);
		for (auto& p : cash_flows)
		{
			p = coupon_value;
		}
		return cash_flows;
	}

	template<typename T>
	template<typename S>
	T Bond<T>::compute_yield(const T& i_price, const S& solver, const DF_type& df_type) const
	{
		assert(solver.ndv == 1);
		auto f = [&,use_penalty_method = solver.use_penalty_method](const auto& solution) { return fitness_irr(solution, i_price, nominal_value, cash_flows, time_periods, df_type, use_penalty_method); };
		auto c = [&,constraints_type = solver.constraints_type](const auto& solution) { return constraints_irr(solution, constraints_type); };
		auto res = solve(f, c, solver, "YTM");
		T yield = res[0];
		return yield;
	}

	template<typename T>
	template<typename S>
	T Bond<T>::compute_yield(const T& i_price, const S& solver, const DF_type& df_type, const std::string& bonds_identifier) const
	{
		assert(solver.ndv == 1);
		auto f = [&, use_penalty_method = solver.use_penalty_method](const auto& solution) { return fitness_irr(solution, i_price, nominal_value, cash_flows, time_periods, df_type, use_penalty_method); };
		auto c = [&, constraints_type = solver.constraints_type](const auto& solution) { return constraints_irr(solution, constraints_type); };
		std::string problem = "YTM";
		auto res = solve(f, c, solver, problem.append(bonds_identifier));
		T yield = res[0];
		return yield;
	}

	template<typename T>
	T Bond<T>::compute_macaulay_duration(const DF_type& df_type) const
	{
		assert(yield > 0 && yield < 1);
		assert(cash_flows.size() > 0);
		assert(nominal_value > 0);
		assert(frequency > 0);
		//! Discount factor 
		T discount_factor = 0.0;
		//! Prest cash flows
		T denominator = 0.0;
		//! Present value
		T pv = 0.0;
		//! Macaulay duration
		T numerator = 0.0;
		for (size_t i = 0; i < time_periods.size(); ++i)
		{
			discount_factor = compute_discount_factor(yield, time_periods[i], df_type);
			pv = coupon_value * discount_factor;
			numerator = numerator + pv * time_periods[i];
			denominator = denominator + pv;
		}
		pv = nominal_value * discount_factor;
		numerator = numerator + pv * time_periods.back();
		denominator = denominator + pv;
		T duration = numerator / denominator;
		return duration;
	}
}