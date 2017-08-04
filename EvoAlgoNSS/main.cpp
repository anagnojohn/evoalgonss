// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "svensson.h"
#include "bond.h"
#include "irr.h"
#include "bondpricing.h"
#include "dependencies.h"
#include "benchmarks.h"

template<typename T>
std::vector<Bond<T>> read_bonds_from_file(const std::string & filename)
{
	std::vector<Bond<T>> bonds;
	std::ifstream input(filename);
	for (std::string line; getline(input, line); )
	{
		T coupon_percentage;
		T price;
		T nominal_value;
		T frequency;
		std::string settlement_date;
		std::string maturity_date;
		std::istringstream stream(line);
		stream >> coupon_percentage >> price >> nominal_value >> frequency >> settlement_date
			>> maturity_date;
		const Bond<T> bond { coupon_percentage, price, nominal_value, frequency, settlement_date, maturity_date };
		bonds.push_back(bond);
	}
	return bonds;
}

int main()
{
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	auto ir_vec = read_ir_from_file<double>("interest_rate_data_periods.txt");
	auto bonds = read_bonds_from_file<double>("bond_data.txt");
	// IRR solvers
	DEstruct<double> de_irr{ 0.6, 1,{ 0.5 },{ 0.7 }, 200, 0.0000001, 200 };
	GAstruct<double> ga_irr { 0.4, 0.35, { 0.5 }, { 0.7 }, 200, 0.0000001, 200 };
	PSOstruct_inertia<double> pso_inertia_irr{ 1.0, 1.0, 5, 0.9, { 0.5 }, { 0.7 }, 200, 0.0000001, 200 };
	PSOstruct_clamping<double> pso_clamping_irr{ 1.0, 1.0, 5, 2.0, { 1000 }, { 0.5 }, { 0.7 }, 200, 0.0000001, 200 };
	// Pricing solvers
	auto decision_variables = set_init_nelson_param(bonds, de_irr);
	DEstruct<double> de_pricing{ 0.6, 1, decision_variables, stdev, 60, 0.01, 400 };
	DEstruct<double> de_fitting{ 0.6, 1, decision_variables, stdev, 60, 0.0001, 400 };
	GAstruct<double> ga_pricing { 0.4, 0.35, decision_variables, stdev_ga, 800, 0.00001, 400 };
	PSOstruct_clamping<double> pso_clamping_pricing { 1.0, 1.0, 5, 2.0, { 1000, 1000, 1000, 1000, 1000, 1000}, decision_variables, stdev, 200, 0.001, 200 };
	PSOstruct_inertia<double> pso_inertia_pricing{ 1.0, 1.0, 5, 0.9, decision_variables, stdev, 200, 0.001, 500 };
	//std::vector <double> cash_flows{ 10.0, 10.0, 10.0, 10.0 };
	//std::cout << macaulay_duration(0.04, cash_flows, 100.0, 2.0);
	//std::cout << macaulay_duration2(0.04, cash_flows, 100.0, 2.0);
	// Call benchmark functions
	benchmarkbondpricing(bonds, de_pricing);
	benchmarkbondpricing_yields(bonds, de_pricing);
	//benchmarkyieldcurvefitting(ir_vec, de_fitting);
	//benchmarkbondpricing(bonds, ga_pricing);
	benchmarkyieldcurvefitting(ir_vec, ga_pricing);
	//benchmarkbondpricing(bonds, pso_clamping_pricing);
	//benchmarkbondpricing(bonds, pso_inertia_pricing);
    return 0;
}