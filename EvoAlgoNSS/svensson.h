#pragma once

#include "bond.h"


//! Nelson-Siegel-Svensson (NSS) model
namespace nss
{
	//! Constraints function for the NSS model
	template<typename T>
	bool constraints_svensson(const std::vector<T>& solution)
	{
		const T& b0 = solution[0];
		const T& b1 = solution[1];
		const T& b2 = solution[2];
		const T& b3 = solution[3];
		const T& tau1 = solution[4];
		const T& tau2 = solution[5];
		if (b0 > 0 && b0 + b1 > 0 && tau1 > 0 && tau2 > 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	//! Spot interest rate at term m using the NSS model
	template<typename T>
	T svensson(const std::vector<T>& solution, const T& m)
	{
		const T& b0 = solution[0];
		const T& b1 = solution[1];
		const T& b2 = solution[2];
		const T& b3 = solution[3];
		const T& tau1 = solution[4];
		const T& tau2 = solution[5];
		if (m == 0)
		{
			return b0 + b1;
		}
		else
		{
			T result = b0 + (b1 + b2) * (tau1 / m) * (1 - std::exp(-m / tau1))
				- b2 * std::exp(-m / tau1) + b3 * (tau2 / m) * (1 - std::exp(-m / tau2))
				- b3 * std::exp(-m / tau2);
			return result;
		}
	}

	//! Penalty function for NSS
	template<typename T>
	T penalty_svensson(const std::vector<T>& solution)
	{
		const T& b0 = solution[0];
		const T& b1 = solution[1];
		const T& b2 = solution[2];
		const T& b3 = solution[3];
		const T& tau1 = solution[4];
		const T& tau2 = solution[5];
		const T C = 1000;
		T sum = 0;
		if (b0 < 0 || b0 > 15)
		{
			sum = sum + C * std::pow(std::abs(b0), 2);
		}
		if (b0 + b1 < 0)
		{
			sum = sum + C * std::pow(std::abs(b0 + b1), 2);
		}
		if (b1 < -15 || b1 > 30)
		{
			sum = sum + C * std::pow(std::abs(b1), 2);
		}
		if (b2 < -30 || b2 > 30)
		{
			sum = sum + C * std::pow(std::abs(b1), 2);
		}
		if (b3 < -30 || b3 > 30)
		{
			sum = sum + C * std::pow(std::abs(b1), 2);
		}
		if (tau1 < 0 || tau1 > 2.5)
		{
			sum = sum + C * std::pow(std::abs(tau2), 2);
		}
		if (tau2 < 2.5 || tau2 > 5.5)
		{
			sum = sum + C * std::pow(std::abs(tau2), 2);
		}
		return sum;
	}
}