//! ConsoleApplication1.cpp : Defines the entry point for the console application.
//!

#include "stdafx.h"
#include "svensson.h"
#include "bond.h"
#include "bondhelper.h"

int main()
{
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	double irr_tol = 0.0000001;
	//! IRR solvers
	DE<double> de_irr{ 0.6, 1, { 0.05 },{ 0.7 }, 30, irr_tol, 200 };
	GA<double> ga_irr{ 0.4, 0.35, 6.0, { 0.5 }, { 0.7 }, 200, irr_tol, 200 };
	PSO<double> pso_irr{ 2.05, 2.05, 4, 0.729, 1, { 1000000 }, { 0.05 }, { 0.7 }, 16, irr_tol, 500 };
	//! Call benchmark functions
	InterestRate_Helper<double> ir{ read_ir_from_file<double>("interest_rate_data_periods.txt") };
	BondHelper<double> de{ read_bonds_from_file<double>("bond_data.txt") };
	auto decision_variables = de.set_init_nss_params(de_irr);
	//de.set_init_nss_params(ga_irr);
	de.set_init_nss_params(pso_irr);
	//! Pricing solvers
	DE<double> de_pricing{ 0.6, 1, decision_variables, stdev, 60, 0.01, 400 };
	DE<double> de_fitting{ 0.6, 1, decision_variables, stdev, 60, 0.0001, 400 };
	GA<double> ga_pricing{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 800, 0.00001, 400 };
	//Local_Best_PSO<double> pso_pricing{ 2.05, 2.05, 6, 0.729, 1.0, { 1000000, 1000000, 1000000, 1000000, 1000000, 1000000 }, decision_variables, stdev, 24, 0.001, 2000};
	//de.bondpricing_prices(de_pricing);
	//de.bondpricing_yields(de_pricing);
	//yieldcurve_fitting(ir_vec, de_fitting);
	//de.bondpricing_yields(ga_pricing);
	//!yieldcurve_fitting(ir_vec, ga_pricing);
	//!de.bondpricing_prices(pso_pricing);
	//!de.bondpricing_yields(pso_pricing);
    return 0;
}