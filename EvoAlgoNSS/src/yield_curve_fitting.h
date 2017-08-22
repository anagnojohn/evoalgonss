#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "svensson.h"

using namespace nss;

//! Yield Curve Fitting namespace
namespace yft {
	/*! \struct Interest_Rate
	*  \brief Structure for interest rates
	*/
	template<typename T>
	struct Interest_Rate
	{
		/** \fn Interest_Rate(const T& i_period, const T& i_rate)
		*  \brief Constructor
		*  \param i_period The period the zero rate was recorded
		*  \param i_rate Zero rate (spot interest rate)
		*  \return An Interest_Rate object
		*/
		Interest_Rate(const T& i_period, const T& i_rate) :
			period{ i_period },
			rate{ i_rate } {};
		/** \brief Period is the time when the rate was recorder */
		const T period;
		/** \brief Interest rate */
		const T rate;
	};

	/** \fn read_ir_from_file(const std::string & filename)
	*  \brief Reads the interest rates and periods from file and constructs a vector of interest rate structs
	*  \param filename The name of the input file as an std::string
	*  \return A vector of Interest_Rate objects
	*/
	template<typename T>
	std::vector<Interest_Rate<T>> read_ir_from_file(const std::string & filename)
	{
		std::vector<Interest_Rate<T>> ir_vec;
		std::ifstream input(filename);
		for (std::string line; getline(input, line); )
		{
			T period;
			T rate;
			std::istringstream stream(line);
			stream >> period >> rate;
			const Interest_Rate<T> ir{ period, rate };
			ir_vec.push_back(ir);
		}
		return ir_vec;
	}

	/*! \class Interest_Rate_Helper
	*  \brief  A class for the yield-curve-fitting problem
	*/
	template<typename T>
	class Interest_Rate_Helper
	{
	public:
		/** \fn Interest_Rate_Helper(const std::vector<Interest_Rate<T>>& i_ir_vec)
		*  \brief Constructor
		*  \param i_ir_vec A vector of Interest_Rate<T> objects
		*  \return An Interest_Rate_Helper object
		*/
		Interest_Rate_Helper(const std::vector<Interest_Rate<T>>& i_ir_vec) :
			ir_vec{ i_ir_vec } {};
		/** \fn yieldcurve_fitting(const S& solver)
		*  \brief Yield Curve Fitting using interest rates and recorded periods
		*  \param solver The parameter structure of the solver that is going to be used for yield curve fitting
		*  \return void
		*/
		template<typename S>
		void yieldcurve_fitting(const S& solver)
		{
			assert(solver.ndv == 6);
			auto f = [&, use_penalty_method = solver.use_penalty_method](const auto& solution) { return fitness_yield_curve_fitting(solution, use_penalty_method); };
			auto c = [&, constraints_type = solver.constraints_type](const auto& solution) { return constraints_svensson(solution, constraints_type); };
			std::cout << "Yield Curve fitting." << "\n";
			auto res = solve(f, c, solver, "YFT");
			T error = 0;
			for (const auto& p : ir_vec)
			{
				error = error + std::pow(svensson(res, p.period) - p.rate, 2);
				//std::cout << "Estimated interest rates: " << svensson(res, p.period) << " Actual interest rates: " << p.rate << "\n";
			}
			std::cout << "Zero-rate Mean Squared Error: " << error / ir_vec.size() << "\n";
		};
	private:
		/** \brief Vector of interest rates */
		std::vector<Interest_Rate<T>> ir_vec;
		/** \fn fitness_yield_curve_fitting(const std::vector<T>& solution, const bool& use_penalty_method)
		*  \brief This is the fitness function for yield-curve fitting using Interest Rates
		*  \param solution NSS parameters candindate solution
		*  \param use_penalty_method Whether to use the penalty method defined for NSS or not
		*  \return The fitness cost of NSS for yield curve fitting
		*/
		T fitness_yield_curve_fitting(const std::vector<T>& solution, const bool& use_penalty_method)
		{
			//! The sum of squares of errors betwwen the actual rates and the rates computed by svensson are used
			T sum_of_squares = 0;
			for (size_t i = 0; i < ir_vec.size(); ++i)
			{
				T estimate = svensson(solution, ir_vec[i].period);
				sum_of_squares = sum_of_squares + std::pow(ir_vec[i].rate - estimate, 2);
			}
			if (use_penalty_method)
			{
				return sum_of_squares + penalty_svensson(solution);
			}
			else
			{
				return sum_of_squares;
			}
		};
	};
}