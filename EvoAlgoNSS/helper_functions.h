#pragma once

//! Overload the operator << for printing vectors
template <typename T>
std::ostream& operator<< (std::ostream& stream, const std::vector<T>& vector) {
	if (!vector.empty()) {
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
	std::uniform_real_distribution<T> i_distribution(0.0, 1.0);
	distribution = i_distribution;
	individuals = init_individuals(decision_variables, npop, stdev);
	min_cost = decision_variables;
	iter = 0;
	population_constraints_checker(decision_variables, stdev, c);
	find_min_cost(f);
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

//! Another version of Macaulay Duration
template<typename T>
T macaulay_duration2(const T& yield, const std::vector<T>& cash_flows, const T& nominal_value, const T& frequency)
{
	assert(yield > 0 && yield < 1);
	assert(cash_flows.size() > 0);
	assert(nominal_value > 0);
	assert(frequency > 0);
	T discount_factor = 0.0;
	T pv_cash_flow = 0.0;
	T duration = 0.0;
	T coupon_value = cash_flows[0];
	T pv = 0.0;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = std::exp(-yield * (static_cast<T>(i + 1) / frequency));
		pv_cash_flow = pv_cash_flow + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + nominal_value * discount_factor;
	for (auto i = 0; i < cash_flows.size(); ++i)
	{
		discount_factor = std::exp(-yield * (static_cast<T>(i + 1) / frequency));
		pv = coupon_value * discount_factor;
		duration = duration + (static_cast<T>(i + 1) / frequency) * pv / pv_cash_flow;
	}
	pv = nominal_value * discount_factor;
	duration = duration + (static_cast<T>(cash_flows.size()) / frequency) * pv / pv_cash_flow;
	return duration;
}