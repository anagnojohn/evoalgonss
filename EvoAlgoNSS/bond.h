#pragma once

#include "dependencies.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"
#include "irr.h"

template<typename T>
class Bond
{
public:
	const T coupon_percentage;
	const T price;
	const T nominal_value;
	const size_t frequency;
	boost::gregorian::date settlement_date;
	boost::gregorian::date maturity_date;
	Bond(const T& i_coupon_percentage, const T& i_price, const T& i_nominal_value, const size_t& i_frequency, const std::string& i_settlement_date, const std::string& i_maturity_date)
		: coupon_percentage(i_coupon_percentage), price(i_price), nominal_value(i_nominal_value), frequency(i_frequency),
		settlement_date(boost::gregorian::from_simple_string(i_settlement_date)),
		maturity_date(boost::gregorian::from_simple_string(i_maturity_date))
	{
		if (settlement_date > maturity_date)
		{
			std::swap(settlement_date, maturity_date);
		}
		int time_periods = maturity_date.year() - settlement_date.year();
		if (maturity_date.month() < settlement_date.month() || (maturity_date.month() == settlement_date.month() && maturity_date.day() < settlement_date.day()))
		{
			time_periods -= 1; 
		}
		const size_t num_time_periods = static_cast<size_t>(time_periods) * frequency;
		cash_flows.resize(num_time_periods);
		for (auto& p : cash_flows)
		{
			p = (coupon_percentage * nominal_value) / frequency;
		}
		coupon_value = coupon_percentage * nominal_value;
		duration = macaulay_duration();
	}
	T yield;
	T duration;
	T coupon_value;
	std::vector<T> cash_flows;
	T macaulay_duration();
private:
	
};

template<typename T>
T Bond<T>::macaulay_duration()
{
	T discount_factor;
	T pv_cash_flow = 0.0;
	T current_bond_price = 0.0;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = 1 / std::pow(1 + coupon_percentage, i + 1);
		pv_cash_flow = pv_cash_flow + (i + 1) * cash_flows[i] * discount_factor;
		current_bond_price = current_bond_price + cash_flows[i] * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + cash_flows.size() * nominal_value * discount_factor;
	current_bond_price = current_bond_price + nominal_value * discount_factor;
	return pv_cash_flow / current_bond_price;
}

template<typename T, typename S>
T getyield(Bond<T> b, S& solver)
{
	auto f = [&](T r) { return irr(r, b.nominal_value, b.cash_flows, b.frequency);};
	return solve(f, b.price, solver)[0];
}

