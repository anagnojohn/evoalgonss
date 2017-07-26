#pragma once

#include <vector>
#include "dependencies.h"
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/time.hpp"

template<typename T>
std::vector<T> compute_cash_flows(T coupon_percentage, T nominal_value, size_t frequency, boost::gregorian::date settlement_date, boost::gregorian::date maturity_date)
{
	//assert(settlement_date < maturity_date);
	assert(coupon_percentage > 0 && coupon_percentage < 1);
	assert(nominal_value > 0);
	assert(frequency > 0);
	size_t number_of_days_coupon = 360 / frequency;
	size_t time_periods = static_cast<size_t>(maturity_date.year() - settlement_date.year());
	if (maturity_date.month() < settlement_date.month() || (maturity_date.month() == settlement_date.month() && maturity_date.day() < settlement_date.day()))
	{
		time_periods -= 1;
	}
	const size_t num_time_periods = static_cast<size_t>(time_periods) * frequency;
	std::vector<T> cash_flows;
	cash_flows.resize(num_time_periods);
	for (auto& p : cash_flows)
	{
		p = (coupon_percentage * nominal_value) / frequency;
	}
	return cash_flows;
}

template<typename T>
class Bond
{
public:

	// In the case price is given but the yield and the macaulay duration are not given
	Bond(const T& i_coupon_percentage, const T& i_price, const T& i_nominal_value, const size_t& i_frequency, const std::string& i_settlement_date, const std::string& i_maturity_date)
		: coupon_percentage{ i_coupon_percentage }, price{ i_price }, nominal_value{ i_nominal_value }, frequency{ i_frequency },
		settlement_date{ boost::gregorian::from_simple_string(i_settlement_date) },
		maturity_date{ boost::gregorian::from_simple_string(i_maturity_date) },
		coupon_value{ coupon_percentage * nominal_value }
	{
		assert(price > 0);
		cash_flows = compute_cash_flows(coupon_percentage, nominal_value, frequency, settlement_date, maturity_date);
	}
	// In case the price is not given but yield and macaulay duration are given
	Bond(const T& i_coupon_percentage, const T& i_yield, const T& i_duration, const T& i_nominal_value, const size_t& i_frequency, const std::string& i_settlement_date, const std::string& i_maturity_date)
		: coupon_percentage{ i_coupon_percentage }, yield{ i_yield }, duration{ i_duration }, nominal_value {i_nominal_value }, frequency{ i_frequency },
		settlement_date{ boost::gregorian::from_simple_string(i_settlement_date) },
		maturity_date{ boost::gregorian::from_simple_string(i_maturity_date) },
		coupon_value { coupon_percentage * nominal_value }
	{
		assert(yield > 0 && yield < 1);
		assert(duration > 0);
		cash_flows = compute_cash_flows(coupon_percentage, nominal_value, frequency, settlement_date, maturity_date);
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
	const size_t frequency;
	const T coupon_value;
	T yield;
	T duration;
	std::vector<T> cash_flows;
	T get_duration() const { return duration; };
private:
	const boost::gregorian::date settlement_date;
	const boost::gregorian::date maturity_date;
};

template<typename T>
T macaulay_duration(const T& yield, const std::vector<T>& cash_flows, const T& nominal_value, const size_t& frequency)
{
	assert(yield > 0 && yield < 1);
	assert(cash_flows.size() > 0);
	assert(nominal_value > 0);
	assert(frequency > 0);
	T discount_factor = 0.0;
	T pv_cash_flow = 0.0;
	T current_bond_price = 0.0;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = 1 / std::pow(1 + (yield / static_cast<T>(frequency)), static_cast<T>(i + 1));
		T coupon_value = cash_flows[i] / static_cast<T>(frequency);
		pv_cash_flow = pv_cash_flow + (i + 1) * coupon_value * discount_factor;
		current_bond_price = current_bond_price + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + cash_flows.size() * nominal_value * discount_factor;
	current_bond_price = current_bond_price + nominal_value * discount_factor;
	return pv_cash_flow / current_bond_price;
}

