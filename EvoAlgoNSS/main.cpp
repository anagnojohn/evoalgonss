// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "svensson.h"
#include "bond.h"
#include "benchmarks.h"

int main()
{
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	auto ir_vec = read_ir_from_file<double>("interest_rate_data_periods.txt");
	auto bonds = read_bonds_from_file<double>("bond_data.txt");
	// IRR solvers
	DEstruct<double> de_irr{ 0.6, 1,{ 0.5 },{ 0.7 }, 200, 0.0000001, 200 };
	GAstruct<double> ga_irr { 0.4, 0.35, 6.0, { 0.5 }, { 0.7 }, 200, 0.0000001, 200 };
	PSOstruct_inertia<double> pso_inertia_irr{ 0.8, 0.8, 10, 0.9, { 0.5 }, { 0.7 }, 100, 0.0000001, 200 };
	PSOstruct_clamping<double> pso_clamping_irr{ 1.0, 1.0, 5, 2.0, { 1000 }, { 0.5 }, { 0.7 }, 200, 0.0000001, 200 };
	// Call benchmark functions
	BondHelper<double> de { bonds };
	auto decision_variables = de.set_init_nss_params(de_irr);
	//auto decision_variables = de.set_init_nss_params(ga_irr);
	de.set_init_nss_params(pso_inertia_irr);
	// Pricing solvers
	DEstruct<double> de_pricing{ 0.6, 1, decision_variables, stdev, 60, 0.01, 400 };
	DEstruct<double> de_fitting{ 0.6, 1, decision_variables, stdev, 60, 0.0001, 400 };
	GAstruct<double> ga_pricing{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 800, 0.00001, 400 };
	PSOstruct_clamping<double> pso_clamping_pricing{ 1.0, 1.0, 5, 2.0,{ 1000, 1000, 1000, 1000, 1000, 1000 }, decision_variables, stdev, 200, 0.001, 200 };
	PSOstruct_inertia<double> pso_inertia_pricing{ 5, 5, 5, 0.9, decision_variables, stdev, 200, 0.001, 500 };
	de.bondpricing_prices(de_pricing);
	de.bondpricing_yields(de_pricing);
	yieldcurve_fitting(ir_vec, de_fitting);
	de.bondpricing_yields(ga_pricing);
	yieldcurve_fitting(ir_vec, ga_pricing);
	de.bondpricing_prices(pso_clamping_pricing);
	de.bondpricing_yields(pso_inertia_pricing);
    return 0;
}