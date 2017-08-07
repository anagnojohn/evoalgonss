#pragma once

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	if (!v.empty()) {
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	return out;
}

//std::vector<T> diff(npop);
//std::vector<T> sum(ndv);
//std::vector<T> mean(ndv);
//std::vector<T> stdev(ndv);
//std::vector<T> epsilon(ndv);
/*
for (auto& p : sum)
{
p = 0.0;
}
for (auto& p : individuals)
{
for (auto j = 0; j < ndv; ++j)
{
sum[j] = sum[j] + p[j];
}
}
for (auto j = 0; j < ndv; ++j)
{
mean[j] = sum[j] / (npop);
}
for (auto j = 0; j < ndv; ++j)
{
for (auto i = 0; i < npop; ++i)
{
diff[i] = individuals[i][j] - mean[j];
}
T sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
stdev[j] = std::sqrt(sq_sum / npop);
}
*/

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