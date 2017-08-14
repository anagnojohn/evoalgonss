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



/*
template<typename T>
template<typename F, typename C>
bool initial_solution_checker(F f, C c)
{
	if (tol > std::abs(fitness_cost))
	{
		T timer = 0;
		display_results();
		std::cout << "Elapsed time in seconds: " << timer << "\n";
		return true;
	}
	else
	{
		return false;
	}
}

template<typename T, typename S>
void Solver_base<T, S>::display_results()
{
	std::cout << "Optimum solution: " << min_cost << " Fitness Value: " << fitness_cost << "\n";
	std::cout << "Population: " << individuals.size() << " Solved at iteration: " << iter << "\n";

}
*/