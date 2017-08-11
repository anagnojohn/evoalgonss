//! ConsoleApplication1.cpp : Defines the entry point for the console application.
//!

#include "stdafx.h"
#include "svensson.h"
#include "bond.h"
#include "bondhelper.h"

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

int main()
{
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	double irr_tol = 0.0000001;
	//! Call benchmark functions
	InterestRate_Helper<double> ir{ read_ir_from_file<double>("interest_rate_data_periods.txt") };
	BondHelper<double> de{ read_bonds_from_file<double>("bond_data.txt") };
	//! IRR solvers
	DE<double> de_irr{ 0.6, 1,{ 0.05 },{ 0.7 }, 30, irr_tol, 200 };
	GA<double> ga_irr{ 0.4, 0.35, 6.0,{ 0.5 },{ 0.7 }, 200, irr_tol, 200 };
	PSO<double> pso_irr{ 2.05, 2.05, 4, 0.729, 1,{ 1000000 },{ 0.05 },{ 0.7 }, 16, irr_tol, 500 };
	auto decision_variables = de.set_init_nss_params(de_irr);
	//de.set_init_nss_params(ga_irr);
	de.set_init_nss_params(pso_irr);
	DE<double> de_fitting{ 0.6, 1, decision_variables, stdev, 60, 0.0001, 400 };
	GA<double> ga_fitting{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 800, 0.00001, 400 };
	ir.yieldcurve_fitting(de_fitting);
	//ir.yieldcurve_fitting(ga_fitting);
	//! Pricing solvers
	DE<double> de_pricing{ 0.6, 1, decision_variables, stdev, 60, 0.01, 400 };
	GA<double> ga_pricing{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 800, 0.00001, 400 };
	PSO<double> pso_pricing{ 2.05, 2.05, 6, 0.729, 1.0, { 1000000, 1000000, 1000000, 1000000, 1000000, 1000000 }, decision_variables, stdev, 24, 0.001, 2000};
	de.bondpricing_prices(de_pricing);
	de.bondpricing_yields(de_pricing);
	de.bondpricing_yields(ga_pricing);
	de.bondpricing_prices(pso_pricing);
	de.bondpricing_yields(pso_pricing);
    return 0;
}