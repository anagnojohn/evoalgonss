#pragma once

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

//! Calculates the bond price
template<typename T>
T find_bond_price(const T& ytm, const T& coupon_value, const T& nominal_value, const std::vector<T>& time_periods)
{
	const size_t& num_time_periods = time_periods.size();
	T sum = 0.0;
	T pv_cash_flow = 0.0;
	T discount_factor = 0.0;
	for (auto i = 0; i < num_time_periods; ++i)
	{
		discount_factor = std::exp(-ytm * time_periods[i]);
		pv_cash_flow = pv_cash_flow + coupon_value * discount_factor;
	}
	pv_cash_flow = pv_cash_flow + nominal_value * discount_factor;
	return pv_cash_flow;
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

//! This method sets the nss initial svensson parameters using test data
template<typename T>
std::vector<T> set_init_nss_params(std::vector<Bond<T>>& bonds)
{
	read_bonds_from_file<double>("bond_data.txt")
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