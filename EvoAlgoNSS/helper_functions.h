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