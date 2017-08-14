#pragma once

#include <iostream>
#include <vector>

//! Overload the operator << for printing vectors
template <typename T>
std::ostream& operator<<(std::ostream& stream, const std::vector<T>& vector)
{
	if (!vector.empty())
	{
		stream << '[';
		std::copy(vector.begin(), vector.end(), std::ostream_iterator<T>(stream, ", "));
		stream << "\b\b]";
	}
	return stream;
}