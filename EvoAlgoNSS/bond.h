#pragma once

#include <vector>
#include "irr.h"
#include "dependencies.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/time.hpp"

template<typename T>
class BondHelper;

// Calculate the cash flows of the bond
template<typename T>
std::vector<T> compute_cash_flows(const T& coupon_value, const T& frequency, boost::gregorian::date settlement_date, boost::gregorian::date maturity_date)
{
	//assert(settlement_date < maturity_date);
	const T number_of_days_coupon = 365 / frequency;
	const auto days_difference = maturity_date.day_number() - settlement_date.day_number();
	const auto time_periods = static_cast<T>(days_difference) / number_of_days_coupon;
	const size_t num_time_periods = static_cast<size_t>(std::ceil(time_periods));
	std::vector<T> cash_flows(num_time_periods);
	for (auto& p : cash_flows)
	{
		p = coupon_value;
	}
	return cash_flows;
}

//! Bond Class
template<typename T>
class Bond
{
	//! Friend class BondHelper
	friend class BondHelper<T>;
public:
	//! Construct in the case price is given but the yield and the macaulay duration are not given
	Bond(const T& i_coupon_percentage, const T& i_price, const T& i_nominal_value, const T& i_frequency, const std::string& i_settlement_date, const std::string& i_maturity_date)
		: coupon_percentage{ i_coupon_percentage }, price{ i_price }, nominal_value{ i_nominal_value }, frequency{ i_frequency },
		settlement_date{ boost::gregorian::from_simple_string(i_settlement_date) },
		maturity_date{ boost::gregorian::from_simple_string(i_maturity_date) },
		coupon_value{ coupon_percentage * nominal_value / frequency }
	{
		assert(price > 0);
		assert(coupon_percentage > 0 && coupon_percentage < 1);
		assert(nominal_value > 0);
		assert(frequency > 0);
		cash_flows = compute_cash_flows(coupon_value, frequency, settlement_date, maturity_date);
		time_periods.resize(cash_flows.size());
		for (auto i = 0; i < time_periods.size(); ++i)
		{
			time_periods[i] = static_cast<T>(i + 1) / frequency;
		}
	}
	//! Construct in case the price is not given but yield and macaulay duration are given
	Bond(const T& i_coupon_percentage, const T& i_yield, const T& i_duration, const T& i_nominal_value, const T& i_frequency, const std::string& i_settlement_date, const std::string& i_maturity_date)
		: coupon_percentage{ i_coupon_percentage }, yield{ i_yield }, duration{ i_duration }, nominal_value {i_nominal_value }, frequency{ i_frequency },
		settlement_date{ boost::gregorian::from_simple_string(i_settlement_date) },
		maturity_date{ boost::gregorian::from_simple_string(i_maturity_date) },
		coupon_value { coupon_percentage * nominal_value / frequency }
	{
		assert(yield > 0 && yield < 1);
		assert(duration > 0);
		assert(coupon_percentage > 0 && coupon_percentage < 1);
		assert(nominal_value > 0);
		assert(frequency > 0);
		cash_flows = compute_cash_flows(coupon_value, frequency, settlement_date, maturity_date);
		time_periods.resize(cash_flows.size);
		for (auto i = 0; i < time_periods.size(); ++i)
		{
			time_periods[i] = static_cast<T>(i + 1) / frequency;
		}
	}
	//! Construct In the case of yield curve fitting only yield to maturity and macaulay duration is needed
	Bond(const T& i_yield, const T& i_duration)
		: yield{ i_yield }, duration{ i_duration }
	{
		assert(yield > 0 && yield < 1);
		assert(duration > 0);
	}
	//! Bond's annual coupon rate
	const T coupon_percentage;
	//! Bond's price
	const T price;
	//! Bond's face value
	const T nominal_value;
	//! Bond's coupon payment frequency
	const T frequency;
	//! This is the annual coupon divided by the frequency
	const T coupon_value;
	//! Coupon payment periods
	std::vector<T> time_periods;
	//! A vector with all the coupon payments corresponding to time periods
	std::vector<T> cash_flows;
	//! Returns the Macaulday duration of a bond
	T ret_duration() const { return duration; };
	//! Returns the yield-to-maturity of a bond
	T ret_yield() const { return yield; };
	//! Sets the yield and duration if needed
	void set_yield(const T& yield, const T& duration) { yield->this.yield; duration->this.duration; };
	//! Calculates the yield-to-maturity and Macaulay duration using the supplied solver
	template<typename S> void compute_yield(S& solver);
private:
	//! Settlement date of the bond
	const boost::gregorian::date settlement_date;
	//! Maturity date of the bond
	const boost::gregorian::date maturity_date;
	//! Calculates the Macaulay duration of the bond
	T compute_macaulay_duration();
	//! Yield-to-maturity of the bond
	T yield;
	//! Macaulay duration of the bond
	T duration;
};

template<typename T>
template<typename S>
void Bond<T>::compute_yield(S& solver)
{
	assert(solver.ndv == 1);
	auto f = [&](const auto& solution) { return fitness_irr(solution, price, nominal_value, cash_flows, frequency); };
	auto c = [&](const auto& solution) { return constraints_irr(solution); };
	auto res = solve(f, c, solver);
	yield = res[0];
	duration = compute_macaulay_duration();
	std::cout << "Macaulay Duration: " << duration << "\n";
}

// Calculate the macaulay_duration of a bond
template<typename T>
T Bond<T>::compute_macaulay_duration()
{
	assert(yield > 0 && yield < 1);
	assert(cash_flows.size() > 0);
	assert(nominal_value > 0);
	assert(frequency > 0);
	//! Discount factor 
	T discount_factor = 0.0;
	//! Prest cash flows
	T pv_cash_flow = 0.0;
	//! Present value
	T pv = 0.0;
	//! Macaulay duration
	T duration = 0.0;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = 1 / std::pow(1 + (yield / frequency), static_cast<T>(i + 1));
		pv_cash_flow = pv_cash_flow + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + nominal_value * discount_factor;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = 1 / std::pow(1 + (yield / frequency), static_cast<T>(i + 1));
		pv = coupon_value * discount_factor;
		duration = duration + (static_cast<T>(i + 1) / frequency) * pv / pv_cash_flow;
	}
	pv = nominal_value * discount_factor;
	duration = duration + (static_cast<T>(cash_flows.size()) / frequency) * pv / pv_cash_flow;
	return duration;
}

//! Calculates the bond price
template<typename T>
T find_bond_price(const T& ytm, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods)
{
	const size_t& num_time_periods = time_periods.size();
	T sum = 0.0;
	T pv_cash_flow = 0.0;
	T discount_factor = 0.0;
	for (auto i = 0; i < num_time_periods; ++i)
	{
		discount_factor = std::exp(-ytm * time_periods[i]);
		pv_cash_flow = pv_cash_flow + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + nominal_value * discount_factor;
	return pv_cash_flow;
}