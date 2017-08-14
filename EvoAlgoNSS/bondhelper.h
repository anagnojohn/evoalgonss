#pragma once

#include <vector>
#include <tuple>
#include <limits>
#include "bond.h"
#include "svensson.h"
#include "irr.h"
#include "yield_curve_fitting.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"


using namespace ea;
using namespace nss;
namespace bond
{
	//! Reads bonds from file
	template<typename T>
	std::vector<Bond<T>> read_bonds_from_file(const std::string & filename)
	{
		std::vector<Bond<T>> bonds;
		std::ifstream input(filename);
		for (std::string line; getline(input, line); )
		{
			T coupon_percentage;
			T price;
			T nominal_value;
			T frequency;
			std::string settlement_date;
			std::string maturity_date;
			std::istringstream stream(line);
			stream >> coupon_percentage >> price >> nominal_value >> frequency >> settlement_date >> maturity_date;
			const Bond<T> bond{ coupon_percentage, price, nominal_value, frequency, settlement_date, maturity_date };
			bonds.push_back(bond);
		}
		return bonds;
	}

	//! A class for the bond pricing problem
	template<typename T>
	class BondHelper
	{
	public:
		//! Constructor
		BondHelper(const std::vector<Bond<T>>& i_bonds) : bonds{ i_bonds } {};
		//! This method sets the nss initial svensson parameters by computing the bond yields-to-maturity and Macaulay durations
		template<typename S> std::vector<T> set_init_nss_params(S& solver);
		//! This methods solves the bond pricing problem using prices and the supplied solver
		template<typename S> void bondpricing_prices(const S& solver);
		//! This methods solves the bond pricing problem using yield-to-maturities and the supplied solver
		template<typename S> void bondpricing_yields(const S& solver);
		//! This method prints to screen the bond pricing results
		void print_bond_pricing_results(const std::vector<T>& res)
		{
			for (const auto& p : bonds)
			{
				std::cout << "Estimated yield: " << svensson(res, p.duration) << " Actual Yield: " << p.yield << "\n";
			}
			for (const auto& p : bonds)
			{
				std::cout << "Estimated price: " << estimate_bond_pricing(res, p.coupon_value, p.nominal_value, p.time_periods) << " Actual Price: " << p.price << "\n";
			}
		};
	private:
		//! Vector of bonds
		std::vector<Bond<T>> bonds;
		//! Calculates the discount factors
		T discount_factor(T yield, T period);
		//! Returns the bond prices using the estimated spot interest rates computed with svensson
		T estimate_bond_pricing(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods);
		//! This is the fitness function for bond pricing using the bonds' yields-to-maturity
		T fitness_bond_pricing_yields(const std::vector<T>& solution);
		//! This is the fitness function for bond pricing using the bonds' prices
		T fitness_bond_pricing_prices(const std::vector<T>& solution);
	};

	template<typename T>
	T BondHelper<T>::discount_factor(T yield, T period)
	{
		return std::exp(-yield * period);
	}

	template<typename T>
	T BondHelper<T>::estimate_bond_pricing(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods)
	{
		T sum = 0.0;
		//! Call svensson for period
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
	T BondHelper<T>::fitness_bond_pricing_yields(const std::vector<T>& solution)
	{
		//! The sum of squares of errors between the actual bond yield to maturity and the estimated yield to maturity by svensson is used
		T sum_of_squares = 0;
		for (auto i = 0; i < bonds.size(); ++i)
		{
			T estimate = svensson(solution, bonds[i].duration);
			sum_of_squares = sum_of_squares + std::pow((bonds[i].yield - estimate), 2);
		}
		sum_of_squares = sum_of_squares + penalty_svensson(solution);
		return sum_of_squares;
	}

	template<typename T>
	T BondHelper<T>::fitness_bond_pricing_prices(const std::vector<T>& solution)
	{
		//! The sum of squares of errors between the actual bond price and the estimated price from estimate_bond_pricing
		T sum_of_squares = 0.0;
		for (auto i = 0; i < bonds.size(); ++i)
		{
			T estimate = estimate_bond_pricing(solution, bonds[i].coupon_value, bonds[i].nominal_value, bonds[i].time_periods);
			sum_of_squares = sum_of_squares + std::pow((bonds[i].price - estimate) / bonds[i].nominal_value, 2) / std::sqrt(bonds[i].duration);
		}
		sum_of_squares = sum_of_squares + penalty_svensson(solution);
		return sum_of_squares;
	}

	template<typename T>
	template<typename S>
	std::vector<T> BondHelper<T>::set_init_nss_params(S& solver)
	{
		for (auto i = 0; i < bonds.size(); ++i)
		{
			std::cout << "Processing bond: " << i + 1 << "\n";
			bonds[i].compute_yield(solver);
		}
		size_t minimum_index = 0;
		size_t maximum_index = 0;
		for (auto i = 0; i < bonds.size(); ++i)
		{
			if (bonds[i].duration < bonds[minimum_index].duration)
			{
				minimum_index = i;
			}
			else
			{
			}
			if (bonds[i].duration > bonds[maximum_index].duration)
			{
				maximum_index = i;
			}
			else
			{
			}
		}
		double b0 = bonds[minimum_index].yield;
		double b1 = std::abs(b0 - bonds[maximum_index].yield);
		double b2 = 0;
		double b3 = 0;
		T target = (b0 + b0 + b1) / 2;
		size_t index = 0;
		for (auto i = 0; i < bonds.size(); ++i)
		{
			if (std::abs(target - bonds[i].yield) < std::abs(target - bonds[index].yield))
			{
				index = i;
			}
			else
			{
			}
		}
		double tau1 = bonds[index].duration;
		double tau2 = tau1;
		const std::vector<T> decision_variables{ b0, b1, b2, b3, tau1, tau2 };
		return decision_variables;
	}

	template<typename T>
	template<typename S>
	void BondHelper<T>::bondpricing_prices(const S& solver)
	{
		assert(solver.ndv == 6);
		for (const auto& p : bonds)
		{
			assert(p.yield > 0 && p.yield < 1);
			assert(p.duration > 0);
		}
		auto f = [&](const auto& solution) { return fitness_bond_pricing_prices(solution); };
		auto c = [&](const auto& solution) { return constraints_svensson(solution); };
		std::cout << "Solving bond pricing using bond prices..." << "\n";
		auto res = solve(f, c, solver);
		print_bond_pricing_results(res);
	}

	template<typename T>
	template<typename S>
	void BondHelper<T>::bondpricing_yields(const S& solver)
	{
		assert(solver.ndv == 6);
		for (const auto& p : bonds)
		{
			assert(p.yield > 0 && p.yield < 1);
			assert(p.duration > 0);
		}
		auto f = [&](const auto& solution) { return fitness_bond_pricing_yields(solution); };
		auto c = [&](const auto& solution) { return constraints_svensson(solution); };
		std::cout << "Solving bond pricing using bond yields..." << "\n";
		auto res = solve(f, c, solver);
		print_bond_pricing_results(res);
	}

	//! This method sets the nss initial svensson parameters using test data
	template<typename T>
	std::vector<T> set_init_nss_params(std::vector<Bond<T>>& bonds)
	{
		read_bonds_from_file<double>("bond_data.txt");
		bonds[0].yield = 0.054308895;
		bonds[1].yield = 0.090624152;
		bonds[2].yield = 0.030896968;
		bonds[3].yield = 0.006625537;
		bonds[4].yield = 0.07972484;
		bonds[5].yield = 0.03366204;
		bonds[6].yield = 0.039963969;
		bonds[7].yield = 0.070339142;
		bonds[0].duration = 2.944988754;
		bonds[1].duration = 4.966178711;
		bonds[2].duration = 6.279674883;
		bonds[3].duration = 9.474865358;
		bonds[4].duration = 8.416273259;
		bonds[5].duration = 14.93089635;
		bonds[6].duration = 15.38446779;
		bonds[7].duration = 12.22184684;
		double b0 = bonds[0].yield;
		double b1 = bonds[6].yield - b0;
		double b2 = 0;
		double b3 = 0;
		double tau1 = bonds[0].duration;
		double tau2 = tau1;
		const std::vector<T> decision_variables{ b0, b1, b2, b3, tau1, tau2 };
		return decision_variables;
	}
}