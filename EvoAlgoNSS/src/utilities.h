#pragma once

#include <iostream>
#include <vector>
#include <iterator>

//! Utilities namespace
namespace utilities {
	/** \enum DF_type
	*  \brief Enumeration for discount factor types/methods
	*/
	enum class DF_type
	{
		frac, /*!< Use fractional form to calculate the discount factor */
		exp /*!< Use the exponential form to calculate the discount factor */
	};
	/** \enum Constraints_type
	*  \brief Enumeration for types of constraints for the optimisation problems
	*/
	enum class Constraints_type
	{ 
		normal, /*!< Use the normal constraints */
		tight, /*!< Use tighter constraints */
		none /*!< Ignore constraints*/
	};
	/** \fn operator<<(std::ostream& stream, const std::vector<T>& vector)
	*  \brief Overload the operator << for printing vectors
	*  \param stream An out stream
	*  \param vector A vector
	*  \return An out stream of the vector in the form [...]
	*/
	template <typename T>
	std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vector)
	{
		if (!vector.empty())
		{
			stream << "[ ";
			std::copy(vector.begin(), vector.end(), std::ostream_iterator<T>(stream, " "));
			stream << "]";
		}
		return stream;
	}
}