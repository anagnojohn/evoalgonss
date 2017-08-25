/** \file svensson.h
* \author Ioannis Anagnostopoulos
* \brief Functions for the Nelson-Siegel-Svensson model
*/

#pragma once

//! Nelson-Siegel-Svensson (NSS) model namespace
namespace nss
{
	/** \fn constraints_svensson(const std::vector<T>& solution, const Constraints_type& constraints_type)
	*  \brief Constraints function for the NSS model
	*  \param solution NSS parameters candindate solution
	*  \param constraints_type Type of constraints used
	*  \return True if constraints are satisfied, false otherwise
	*/
	template<typename T>
	bool constraints_svensson(const std::vector<T>& solution, const Constraints_type& constraints_type)
	{
		switch (constraints_type)
		{
		case(Constraints_type::normal):
		{
			const T& b0 = solution[0];
			const T& b1 = solution[1];
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
		case(Constraints_type::tight):
		{
			const T& b0 = solution[0];
			const T& b1 = solution[1];
			const T& b2 = solution[2];
			const T& b3 = solution[3];
			const T& tau1 = solution[4];
			const T& tau2 = solution[5];
			if ((b0 > 0 && b0 < 15)
				&& (b1 > -15 && b1 < 30)
				&& (b2 > -30 && b2 < 30)
				&& (b3 > -30 && b3 < 30)
				&& (tau1 > 0 && tau1 < 2.5)
				&& (tau2 > 2.5 && tau2 < 5.5))
			{
				return true;
			}
			else
			{
				return false;
			}
		}
		case(Constraints_type::none): return true;
		default: std::abort();
		}
	}

	/** \fn svensson(const std::vector<T>& solution, const T& m)
	*  \brief Spot interest rate at term m using the NSS model
	*  \param solution Candidate solution for the parameters of NSS
	*  \param m The term at which the spot interest rate is recorded
	*  \return The spot interest rate at term m
	*/
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
			T result = b0;
			result = result + b1 * ((1 - (std::exp(-m / tau1))) / (m / tau1));
			result = result + b2 * (((1 - (std::exp(-m / tau1))) / (m / tau1)) - std::exp(-m / tau1));
			result = result + b3 * (((1 - (std::exp(-m / tau2))) / (m / tau2)) - std::exp(-m / tau2));
			return result;
		}
	}

	/** \fn penalty_svensson(const std::vector<T>& solution)
	*  \brief Penalty function for NSS
	*  \param solution Candidate solution for the parameters of NSS
	*  \return A penalty value, if constraints are not satisfied
	*/
	template<typename T>
	T penalty_svensson(const std::vector<T>& solution)
	{
		T sum = 0;
		const T& b0 = solution[0];
		const T& b1 = solution[1];
		const T& b2 = solution[2];
		const T& b3 = solution[3];
		const T& tau1 = solution[4];
		const T& tau2 = solution[5];
		const T C = 100000;
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