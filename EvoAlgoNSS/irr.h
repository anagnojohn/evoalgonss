#pragma once

template<typename T>
T irr(T r, T future_value, std::vector<T> cash_flows)
{
	const size_t& num_time_periods = cash_flows.size();
	T sum = 0.0;
	for (auto i = 0; i < num_time_periods)
	{
		sum = sum + cash_flows[i];
	}
	return cash_flows / pow((1 + r), num_time_periods) + future_value / pow((1 + r), num_time_periods);
}