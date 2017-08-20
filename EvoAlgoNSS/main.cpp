//! ConsoleApplication1.cpp : Defines the entry point for the console application.
//!

#include "stdafx.h"
#include "options.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"
#include "bondhelper.h"

int main()
{
	using namespace bond;
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	double irr_tol = 0.0000001;
	//! Call benchmark functions
	InterestRate_Helper<double> ir{ read_ir_from_file<double>("interest_rate_data_periods.txt") };
	BondHelper<double> de{ read_bonds_from_file<double>("bond_data_3.txt") };
	//! IRR solvers
	DE<double> de_irr{ 0.6, 1,{ 0.05 },{ 0.7 }, 30, irr_tol, 200, false };
	GA<double> ga_irr{ 0.4, 0.35, 6.0, { 0.05 },{ 0.5 }, 10, irr_tol, 30, false, Constraints_type::normal, Strategy::keep_same, false};
	PSO<double> pso_irr{ 2.05, 2.05, 4, 0.729, 1,{ 1000000 },{ 0.05 },{ 0.7 }, 16, irr_tol, 500, false};
	auto decision_variables = de.set_init_nss_params(de_irr);
	de.set_init_nss_params(ga_irr);
	de.set_init_nss_params(pso_irr);
	DE<double> de_fitting{ 0.6, 0.99, decision_variables, stdev, 60, 0.0001, 1000};
	GA<double> ga_fitting{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 10, 0.00001, 100, false, Constraints_type::normal, Strategy::keep_same };
	ir.yieldcurve_fitting(de_fitting);
	ir.yieldcurve_fitting(ga_fitting);
	//! Pricing solvers
	DE<double> de_pricing{ 0.6, 0.99, decision_variables, stdev, 40, 0.01, 1000 };
	GA<double> ga_pricing{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 800, 0.00001, 400, false, Constraints_type::normal, Strategy::keep_same };
	PSO<double> pso_pricing{ 2.05, 2.05, 6, 0.729, 1.0, { 1.5, 4.5, 6, 6, 0.25, 0.3 }, decision_variables, stdev, 24, 0.001, 2000};
	//ir.yieldcurve_fitting(pso_pricing);
	de.bond_pricing(de_pricing, de_irr, Bond_pricing_type::bpp);
	de.bond_pricing(ga_pricing, ga_irr, Bond_pricing_type::bpp);
	de.bond_pricing(pso_pricing, pso_irr, Bond_pricing_type::bpp);
    return 0;
}