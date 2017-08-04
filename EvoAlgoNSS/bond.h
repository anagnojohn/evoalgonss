#pragma once

#include <vector>
#include "dependencies.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/time.hpp"

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

template<typename T>
class Bond
{
public:

	// In the case price is given but the yield and the macaulay duration are not given
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
	// In case the price is not given but yield and macaulay duration are given
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
	// In the case of yield curve fitting only yield to maturity and macaulay duration is needed
	Bond(const T& i_yield, const T& i_duration)
		: yield{ i_yield }, duration{ i_duration }
	{
		assert(yield > 0 && yield < 1);
		assert(duration > 0);
	}
	const T coupon_percentage;
	const T price;
	const T nominal_value;
	const T frequency;
	const T coupon_value;
	std::vector<T> time_periods;
	T yield;
	T duration;
	std::vector<T> cash_flows;
	T get_duration() const { return duration; };
private:
	const boost::gregorian::date settlement_date;
	const boost::gregorian::date maturity_date;
};

template<typename T>
T macaulay_duration(const T& yield, const std::vector<T>& cash_flows, const T& nominal_value, const T& frequency)
{
	assert(yield > 0 && yield < 1);
	assert(cash_flows.size() > 0);
	assert(nominal_value > 0);
	assert(frequency > 0);
	T discount_factor = 0.0;
	T pv_cash_flow = 0.0;
	//T current_bond_price = 0.0;
	T duration = 0.0;
	T coupon_value = cash_flows[0];
	T pv = 0.0;
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

template<typename T>
T macaulay_duration2(const T& yield, const std::vector<T>& cash_flows, const T& nominal_value, const T& frequency)
{
	assert(yield > 0 && yield < 1);
	assert(cash_flows.size() > 0);
	assert(nominal_value > 0);
	assert(frequency > 0);
	T discount_factor = 0.0;
	T pv_cash_flow = 0.0;
	T duration = 0.0;
	T coupon_value = cash_flows[0];
	T pv = 0.0;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = std::exp(-yield * (static_cast<T>(i + 1) / frequency));
		pv_cash_flow = pv_cash_flow + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + nominal_value * discount_factor;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = std::exp(-yield * (static_cast<T>(i + 1) / frequency));
		pv = coupon_value * discount_factor;
		duration = duration + (static_cast<T>(i + 1) / frequency) * pv / pv_cash_flow;
	}
	pv = nominal_value * discount_factor;
	duration = duration + (static_cast<T>(cash_flows.size()) / frequency) * pv / pv_cash_flow;
	return duration;
}

