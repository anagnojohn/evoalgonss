/** \file main.cpp
* \author Ioannis Anagnostopoulos
* \brief A showcase of the application of the solvers on the Yield Curve Fitting, Internal Rate of Return Estimation and Bond Pricing Problems
* \details Usage is the following:
*
* Step 1:Create a solver structure object {GA, DE, PSOl} with a specific floating-point number type, setting all of its parameters throught its constructor.
*
* Set print_to_output or print_to_display to false if there is no need for displaying the results to terminal or printing them to a file.
*
* Step 2: Either use the common interface solve solve(const F& f, const C& c, const S<T>& solver_struct, const std::string& problem_name) passing
* the objective and constraint functions as lambda functions (anonymous functions) capturing all the required variables, such as bonds
* or use the public interfaces from the Interest_Rate_Helper (Yield Curve Fitting), BondHelper (Internal Rate of Return estimation and bond pricing for a number of bonds)
* and Bond (Internal Rate of Return Estimation and Macaulay Duration Estimation) classes after creating an object instance of those classes using their constructors.
*/

/** \mainpage
* \section intro_sec Introduction
*
* This project implements solvers for three financial optimisation problems: <b>Yield Curve Fitting</b>, <b>Internal Rate of Return Estimation</b> and <b>Bond Pricing</b> problems.
*
* Some general guidelines and remarks are mentioned below.
*
* \section template_guide Template Parameters
*
* A template parameter named <b>T</b> denotes a floating-point number type.
*
* A template parameter named <b>F</b> denotes a lambda function (or function object/functor) type of the objective function
*
* A template parameter named <b>C</b> denotes a lambda function (or function object/functor) type of the constraints function
*
* A template parameter named <b>S</b> denotes a solver structure.
*
* \section compile_guide Compilation
*
* Compilation is possible with the three major compilers, gcc, clang and msvc.
*
* Since some of the features used in the project are from the C++11/C++14 standards, the source has to be compiled with at least <b>-std=c++14</b>.
*
* For gcc and clang, <b>-std=c++17</b> can be used as well, as the source is compatible with the newest C++ ISO standard.
*
* \subsection external_libraries External Libraries
*
* The external libraries used by the project are: <b>Boost</b> (http://www.boost.org/) and <b>date</b> (https://github.com/HowardHinnant/date).
*
* Boost is used for the Beta Distribution implementation and the pi constant.
*
* date is used for date handling of the settlement and maturity dates.
*
* Both can be included only as their header versions.
* 
* \section test_guide Showcase
*
* A showcase of the project is provided under <b>tests/main.cpp</b>.
*
* Be sure that the resulting executable will be run in the same working directory as the data files, otherwise the executable will crash abruptly,
* since exceptions are not implemented yet.
*
* If other data files are used, be sure that their content has the correct format.
*
* For bond data files, the input data have to be of the form: <b>coupon rate (in percentage) price nominal value frequency settlement date (in yyyy-mm-dd form) and maturity date (in yyyy-mm-dd form)</b>.
*
* For example: <b>0.06 101.657 100 2 2016-03-30 2019-06-24</b> has the correct format.
*
* For interest rate data files, the input data have to be of the form: <b>period (as a decimal) zero rate (in percentage)</b>.
*
* For example: <b>0.25 0.079573813</b> has the correct format.
*/

#include "../src/bondhelper.h"
#include "../src/geneticalgo.h"
#include "../src/pso_sub_swarm.h"
#include "../src/differentialevo.h"
#include "../src/lbestpso.h"

int main()
{
	using namespace yft;
	using namespace bond;
	const std::vector<double> stdev { 0.7, 0.7, 0.7, 0.7, 0.7, 0.7 };
	const std::vector<double> stdev_ga{ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5 };
	double irr_tol = 0.00000001;
	double tol = 0.0001;
	double tol_f = 0.001;
	//! Call benchmark functions
	Interest_Rate_Helper<double> ir{ read_ir_from_file<double>("interest_rate_data_periods.txt") };
	BondHelper<double> de{ read_bonds_from_file<double>("bond_data_3.txt"), DF_type::exp };
	//! IRR solvers
	DE<double> de_irr{ 1, 0.6,{ 0.05 },{ 0.7 }, 10, irr_tol, 500, false, Constraints_type::normal, true, true };
	DE<double> de_irr_check{ 1, 0.6,{ 0.05 },{ 0.7 }, 10, irr_tol, 500, false, Constraints_type::normal, false, false };
	GA<double> ga_irr{ 0.4, 0.35, 6.0, { 0.05 },{ 0.5 }, 42, irr_tol, 2000, false, Constraints_type::normal, Strategy::remove, true, true};
	PSOl<double> pso_irr{ 1.49618, 0.9, { 1000000 },{ 0.05 },{ 0.7 }, 22, irr_tol, 3000, false, Constraints_type::normal, true, true};
	auto decision_variables = de.set_init_nss_params(de_irr);
	//auto decision_variables1 = de.set_init_nss_params(pso_irr);
	DE<double> de_pricing{ 1, 0.6, decision_variables, stdev, 60, tol, 500, false, Constraints_type::tight, true, true };
	DE<double> de_fitting{ 1, 0.6, decision_variables, stdev, 60, tol_f, 500, false, Constraints_type::tight, true, true };
	GA<double> ga_pricing{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 250, tol, 2000, false, Constraints_type::tight, Strategy::remove, true, true };
	GA<double> ga_fitting{ 0.4, 0.35, 6.0, decision_variables, stdev_ga, 250, tol_f, 2000, false, Constraints_type::tight, Strategy::remove, true, true };
	PSOl<double> pso_pricing{ 1.49618, 0.9, { 100000, 100000, 100000, 100000, 100000, 100000 }, decision_variables, stdev, 130, tol, 3000, false, Constraints_type::tight, true, true };
	PSOl<double> pso_fitting{ 1.49618, 0.9, { 100000, 100000, 100000, 100000, 100000, 100000 }, decision_variables, stdev, 130, tol_f, 3000, false, Constraints_type::tight, true, true };
	//PSOs<double> pso_pricing{ 2.05, 2.05, 6, 0.9, 1.0,{ 100000, 100000, 100000, 100000, 100000, 100000 }, decision_variables, stdev, 24, tol, 1000, false, Constraints_type::none, true, true };
	for (size_t i = 0; i < 100; ++i)
	{
		de.bond_pricing(ga_pricing, de_irr_check, Bond_pricing_type::bpp);
		ir.yieldcurve_fitting(ga_fitting);
		de.bond_pricing(de_pricing, de_irr_check, Bond_pricing_type::bpp);
		ir.yieldcurve_fitting(de_fitting);
		de.bond_pricing(pso_pricing, de_irr_check, Bond_pricing_type::bpp);
		ir.yieldcurve_fitting(pso_fitting);
		de.set_init_nss_params(de_irr);
		de.set_init_nss_params(pso_irr);
		de.set_init_nss_params(ga_irr);
	}
    return 0;
}