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
	// Enumerator for type of bondpricing, using yields or prices
	enum class Bond_pricing_type { bpp, bpy };

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
		BondHelper(const std::vector<Bond<T>>& i_bonds, const DF_type& i_df_type = DF_type::exp) : bonds{ i_bonds }, df_type{ i_df_type } {};
		//! This method sets the nss initial svensson parameters by computing the bond yields-to-maturity and Macaulay durations
		template<typename S> std::vector<T> set_init_nss_params(const S& solver);
		//! This methods solves the bond pricing problem using prices or yields and the supplied solver
		template<typename S1, typename S2> void bond_pricing(const S1& solver, const S2& solver_irr, const Bond_pricing_type& bond_pricing_type);
		//! This method prints to screen the bond pricing results
		template<typename S> void print_bond_pricing_results(const std::vector<T>& res, const S& solver_irr);
	private:
		//! Vector of bonds
		std::vector<Bond<T>> bonds;
		//! Discount Factor type
		const DF_type df_type;
		//! Returns the bond prices using the estimated spot interest rates computed with svensson
		T estimate_bond_pricing(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods);
		//! This is the fitness function for bond pricing using the bonds' yields-to-maturity
		template<typename S> T fitness_bond_pricing_yields(const std::vector<T>& solution, const S& solver_irr, const bool& use_penalty_method);
		//! This is the fitness function for bond pricing using the bonds' prices
		T fitness_bond_pricing_prices(const std::vector<T>& solution, const bool& use_penalty_method);
	};

	template<typename T>
	template<typename S>
	std::vector<T> BondHelper<T>::set_init_nss_params(const S& solver)
	{
		for (size_t i = 0; i < bonds.size(); ++i)
		{
			std::cout << "Processing bond: " << i + 1 << "\n";
			bonds[i].yield = bonds[i].compute_yield(bonds[i].price, solver, df_type, std::to_string(i+1));
			bonds[i].duration = bonds[i].compute_macaulay_duration(df_type);
			std::cout << "Yield to Maturity: " << bonds[i].yield << "\n";
			std::cout << "Macaulay Duration: " << bonds[i].duration << "\n";
			std::cout << "Estimated Price: " << compute_pv(bonds[i].yield, bonds[i].nominal_value, bonds[i].cash_flows, bonds[i].time_periods, df_type)
				<< " Actual Price: " << bonds[i].price << "\n";
		}
		size_t minimum_index = 0;
		size_t maximum_index = 0;
		for (size_t i = 0; i < bonds.size(); ++i)
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
		for (size_t i = 0; i < bonds.size(); ++i)
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
	T BondHelper<T>::estimate_bond_pricing(const std::vector<T>& solution, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods)
	{
		T sum = 0.0;
		//! Call svensson for period
		for (const auto& t : time_periods)
		{
			sum = sum + coupon_value * compute_discount_factor(svensson(solution, t), t, df_type);
		}		
		T m = time_periods.back();
		sum = sum + nominal_value * compute_discount_factor(svensson(solution, m), m, df_type);
		return sum;
	}

	template<typename T>
	T BondHelper<T>::fitness_bond_pricing_prices(const std::vector<T>& solution, const bool& use_penalty_method)
	{
		//! The sum of squares of errors between the actual bond price and the estimated price from estimate_bond_pricing
		T sum_of_squares = 0.0;
		for (const auto& k : bonds)
		{
			T estimate = estimate_bond_pricing(solution, k.coupon_value, k.nominal_value, k.time_periods);
			sum_of_squares = sum_of_squares + std::pow((k.price/100 - estimate/100), 2) / std::sqrt(k.duration);
		}
		if (use_penalty_method)
		{
			return sum_of_squares + penalty_svensson(solution);
		}
		else
		{
			return sum_of_squares;
		}
	}

	template<typename T>
	template<typename S>
	T BondHelper<T>::fitness_bond_pricing_yields(const std::vector<T>& solution, const S& solver_irr, const bool& use_penalty_method)
	{
		//! The sum of squares of errors between the actual bond yield to maturity and the estimated yield to maturity by svensson is used
		T sum_of_squares = 0;
		for (const auto& k : bonds)
		{
			T estimate_price = estimate_bond_pricing(solution, k.coupon_value, k.nominal_value, k.time_periods);
			T estimate = k.compute_yield(estimate_price, solver_irr, df_type);
			sum_of_squares = sum_of_squares + std::pow(k.yield - estimate, 2);
		}
		if (use_penalty_method)
		{
			return sum_of_squares + penalty_svensson(solution);
		}
		else
		{
			return sum_of_squares;
		}
	}
	
	

	template<typename T>
	template<typename S1, typename S2>
	void BondHelper<T>::bond_pricing(const S1& solver, const S2& solver_irr, const Bond_pricing_type& bond_pricing_type)
	{
		assert(solver.ndv == 6);
		for (const auto& p : bonds)
		{
			assert(p.yield > 0 && p.yield < 1);
			assert(p.duration > 0);
		}
		switch (bond_pricing_type)
		{
		case(Bond_pricing_type::bpp):
		{
			auto f = [&, use_penalty_method = solver.use_penalty_method](const auto& solution) { return fitness_bond_pricing_prices(solution, use_penalty_method); };
			auto c = [&, constraints_type = solver.constraints_type](const auto& solution) { return constraints_svensson(solution, constraints_type); };
			std::cout << "Solving bond pricing using bond prices..." << "\n";
			auto res = solve(f, c, solver, "BPP");
			print_bond_pricing_results(res, solver_irr);
			break;
		}
		case(Bond_pricing_type::bpy):
		{
			auto f = [&, use_penalty_method = solver.use_penalty_method](const auto& solution) { return fitness_bond_pricing_yields(solution, solver_irr, use_penalty_method); };
			auto c = [&, constraints_type = solver.constraints_type](const auto& solution) { return constraints_svensson(solution, constraints_type); };
			std::cout << "Solving bond pricing using bond yields..." << "\n";
			auto res = solve(f, c, solver, "BPY");
			print_bond_pricing_results(res, solver_irr);
		}
		}
	}

	template<typename T>
	template<typename S>
	void BondHelper<T>::print_bond_pricing_results(const std::vector<T>& res, const S& solver_irr)
	{
		T error = 0;
		for (const auto& p : bonds)
		{
			T estimate_price = estimate_bond_pricing(res, p.coupon_value, p.nominal_value, p.time_periods);
			T estimate = p.compute_yield(estimate_price, solver_irr, df_type);
			error = error + std::pow(estimate - p.yield, 2);
			//std::cout << "Estimated yield: " << estimate << " Actual Yield: " << p.yield << "\n";
		}
		std::cout << "Yield Mean Squared Error: " << error / bonds.size() << "\n";
		error = 0;
		for (const auto& p : bonds)
		{
			error = error + std::pow(estimate_bond_pricing(res, p.coupon_value, p.nominal_value, p.time_periods) - p.price, 2);
			//std::cout << "Estimated price: " << estimate_bond_pricing(res, p.coupon_value, p.nominal_value, p.time_periods) << " Actual Price: " << p.price << "\n";
		}
		std::cout << "Price Mean Squared Error: " << error / bonds.size() << "\n";
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