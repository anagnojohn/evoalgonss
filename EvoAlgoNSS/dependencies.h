#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iterator>
//#include <utility>
#include <algorithm>
#include <tuple>
#include <limits>
#include <random>
#include <gsl/gsl_sf_bessel.h>
#include <boost/math/distributions.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/time.hpp"

template<typename T>
void init_epsilon(std::vector <std::vector<T>>& individuals, const std::vector<T>& stdev)
{
	std::random_device generator;
	auto npop = individuals.size();
	auto ndv = individuals[0].size();
	std::vector<T> diff(npop);
	//std::vector<T> sum(ndv);
	//std::vector<T> mean(ndv);
	//std::vector<T> stdev(ndv);
	//std::vector<T> epsilon(ndv);
	T epsilon;
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
	for (auto i = 0; i < npop; ++i)
	{
		for (auto j = 0; j < ndv; ++j)
		{
			std::normal_distribution<> ndistribution(0, stdev[j]);
			epsilon = ndistribution(generator);
			individuals[i][j] = individuals[i][j] + epsilon;
		}
	}
}

template<typename T>
std::vector< std::vector<T>> create_individuals(const size_t& npop, const std::vector<T>& decision_variables)
{
	std::vector< std::vector<T> > individuals(npop, std::vector<T>(decision_variables.size()));
	for (auto i = 0; i < npop; ++i)
	{
		individuals[i] = decision_variables;
	}
	return individuals;
}

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	if (!v.empty()) {
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	return out;
}