//! ConsoleApplication1.cpp : Defines the entry point for the console application.
//!
#include "../src/bondhelper.h"
#include "../src/geneticalgo.h"
#include "../src/local_best_pso.h"
#include "../src/differentialevo.h"

int main()
{
	using namespace yft;
	using namespace bond;
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	double irr_tol = 0.0000001;
	double tol = 0.01;
	//! Call benchmark functions
	Interest_Rate_Helper<double> ir{ read_ir_from_file<double>("interest_rate_data_periods.txt") };
	BondHelper<double> de{ read_bonds_from_file<double>("bond_data_3.txt") };
	//! IRR solvers
	
	DE<double> de_irr{ 0.6, 1.00,{ 0.05 },{ 0.7 }, 20, irr_tol, 400, false, Constraints_type::normal, true, true };
	DE<double> de_irr_check{ 0.6, 1.00,{ 0.05 },{ 0.7 }, 20, irr_tol, 400, false, Constraints_type::normal, false, false };
	GA<double> ga_irr{ 0.4, 0.35, 6.0, { 0.05 },{ 0.5 }, 10, irr_tol, 400, false, Constraints_type::normal, Strategy::keep_same, true, true};
	PSO<double> pso_irr{ 2.05, 2.05, 4, 0.729, 1,{ 1000000 },{ 0.05 },{ 0.7 }, 16, irr_tol, 2000, false, Constraints_type::normal, true, true};
	auto decision_variables = de.set_init_nss_params(de_irr);
	DE<double> de_pricing{ 0.6, 1.00, decision_variables, stdev, 60, tol, 400, false, Constraints_type::normal, true, true };
	GA<double> ga_pricing{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 800, tol, 400, false, Constraints_type::normal, Strategy::keep_same, true, true };
	PSO<double> pso_pricing{ 2.05, 2.05, 6, 0.729, 1.0,{ 100000, 1000000, 1000000, 1000000, 1000000, 1000000}, decision_variables, stdev, 24, tol, 2000, false, Constraints_type::normal, true, true };
	for (auto i = 0; i < 100; ++i)
	{
		//std::cout << "NEW RUN" << "\n" << "\n" << "\n";
		de.set_init_nss_params(ga_irr);
		//de.set_init_nss_params(pso_irr);
		//! Pricing solvers
		//ir.yieldcurve_fitting(de_pricing);
		//ir.yieldcurve_fitting(ga_pricing);
		//ir.yieldcurve_fitting(pso_pricing);
		//de.bond_pricing(de_pricing, de_irr_check, Bond_pricing_type::bpp);
		//de.bond_pricing(ga_pricing, de_irr_check, Bond_pricing_type::bpp);
		//de.bond_pricing(pso_pricing, de_irr_check, Bond_pricing_type::bpp);
	}
    return 0;
}