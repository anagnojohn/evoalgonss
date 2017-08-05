#pragma once

#include "dependencies.h"
#include "svensson.h"
#include "irr.h"
#include "bond.h"

template<typename T>
T discount_factor(T yield, T period)
{
	return std::exp(-yield * period);
}

template<typename T>
T estimate_bond_pricing(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods)
{
	T sum = 0.0;
	for (const auto& p : time_periods)
	{
		sum = coupon_value * discount_factor(svensson(solution, p), p);
	}
	T m = time_periods[time_periods.size() - 1];
	T last_payment = nominal_value * discount_factor(svensson(solution, m), m);
	sum = sum + last_payment;
	return sum;
}

template<typename T>
T fitness_bond_pricing_yields(const std::vector<T>& solution, const std::vector<Bond<T>>& bonds)
{
	T sum_of_squares = 0;
	for (auto i = 0; i < bonds.size(); ++i)
	{
		T estimate = svensson(solution, bonds[i].duration);
		sum_of_squares = sum_of_squares + std::pow((bonds[i].yield - estimate), 2);
	}
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
	const T C = 1000;
	if (b0 < 0 || b0 > 15)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0), 2);
	}
	//if (b0 + b1 < 0)
	//{
	//	sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0 + b1), 2);
	//}
	if (b1 < -15 || b1 > 30)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b1), 2);
	}
	if (b2 < -30 || b2 > 30)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b1), 2);
	}
	if (b3 < -30 || b3 > 30)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b1), 2);
	}
	if (tau1 < 0 || tau1 > 2.5)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	if (tau2 < 2.5 || tau2 > 5.5)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	return sum_of_squares;
}

template<typename T>
T fitness_bond_pricing(const std::vector<T>& solution, const std::vector<Bond<T>>& bonds)
{
	T sum_of_squares = 0.0;
	for (auto i = 0; i < bonds.size(); ++i)
	{
		T estimate = estimate_bond_pricing(solution, bonds[i].coupon_value, bonds[i].nominal_value, bonds[i].time_periods);
		sum_of_squares = sum_of_squares + std::pow((bonds[i].price - estimate) / bonds[i].nominal_value, 2) / std::sqrt(bonds[i].duration);
	}
	const T& b0 = solution[0];
	const T& b1 = solution[1];
	const T& b2 = solution[2];
	const T& b3 = solution[3];
	const T& tau1 = solution[4];
	const T& tau2 = solution[5];
	const T C = 500;
	if (b0 < 0 || b0 > 15)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0), 2);
	}
	//if (b0 + b1 < 0)
	//{
	//	sum_of_squares = sum_of_squares + C * std::pow(std::abs(b0 + b1), 2);
	//}
	if (b1 < -15 || b1 > 30)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b1), 2);
	}
	if (b2 < -30 || b2 > 30)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b1), 2);
	}
	if (b3 < -30 || b3 > 30)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(b1), 2);
	}
	if (tau1 < 0 || tau1 > 2.5)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	if (tau2 < 2.5 || tau2 > 5.5)
	{
		sum_of_squares = sum_of_squares + C * std::pow(std::abs(tau2), 2);
	}
	return sum_of_squares;
}