// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"
#include "svensson.h"
#include "datehandler.h"
#include "bond.h"
#include "irr.h"
#include "bondpricing.h"

int main()
{
	//test();
	double b0 = 4.336 / 100;
	double b1 = 4.991 / 100 - b0;
	double b2 = 0;
	double b3 = 0;
	double b4 = 0;
	double tau1 = 1;
	double tau2 = 2;
	const size_t npop = 500;
	const size_t ndv = 6;
	double maturity = 3;
	const double tol = 0.0001;
	std::vector<std::vector<double>> decision_variables(npop, std::vector<double>(ndv));
	for (auto i = 0; i < decision_variables.size(); ++i)
	{
		decision_variables[i][0] = b0;
		decision_variables[i][1] = b1;
		decision_variables[i][2] = b2;
		decision_variables[i][3] = b3;
		decision_variables[i][4] = tau1;
		decision_variables[i][5] = tau2;
	}

	/*
	std::cout << "Genetic Algorithm:" << '\n';
	GeneticAlgo<double, decltype(f)> object1(decision_variables, tol, 200, 0.4, 0.35, stdev);
	std::cout << "Optimum" << '\n';
	auto b = object1.solve(f, 0.0);
	std::cout << b << " " << f(b) << '\n';
	std::cout << "Particle Swarm:" << '\n';
	LocalBestPSO<double, decltype(f)> object2(decision_variables, tol, 200, 1, 1, 2, 0.9, stdev);
	auto c = object2.solve(f, 0.0);
	std::cout << "Optimum" << '\n';
	std::cout << c << " " << f(c) << '\n';
	*/
	size_t npop_irr = 100;
	std::vector< std::vector<double> > init_rate(npop_irr, std::vector<double>(1));
	for (auto i = 0; i < init_rate.size(); ++i)
	{
		init_rate[i][0] = 0.2;
	}
	std::vector<Bond<double>> bonds;
	Bond<double> bond1(0.06, 101.657, 100, 2, "2016-03-30", "2019-06-24");
	Bond<double> bond2(0.09, 99.675, 100, 2, "2016-03-30", "2022-11-22");
	Bond<double> bond3(0.07, 126.888, 100, 2, "2016-03-30", "2024-01-15");
	Bond<double> bond4(0.04, 136.2, 100, 2, "2016-03-30", "2027-07-09");
	Bond<double> bond5(0.08, 100.22, 100, 2, "2016-03-30", "2030-04-28");
	Bond<double> bond6(0.05, 125.98, 100, 2, "2016-03-30", "2039-03-15");
	Bond<double> bond7(0.06, 132.88, 100, 2, "2016-03-30", "2043-03-15");
	Bond<double> bond8(0.11, 150.33, 100, 2, "2016-03-30", "2048-07-28");
	bonds.push_back(bond1);
	bonds.push_back(bond2);
	bonds.push_back(bond3);
	bonds.push_back(bond4);
	bonds.push_back(bond5);
	bonds.push_back(bond6);
	bonds.push_back(bond7);
	bonds.push_back(bond8);

	std::vector<double> stdev(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		stdev[j] = 0.7;
	}
	size_t ndv_irr = 1;
	std::vector<double> stdev_irr(1);
	for (auto j = 0; j < ndv_irr; ++j)
	{
		stdev_irr[j] = 0.7;
	}
	DifferentialEvo<double> irr_de(init_rate, 0.000000001, 200, 0.5, 0.9, stdev_irr);
	for (auto i = 0; i < bonds.size(); ++i)
	{
		bonds[i].yield = setyield<double, DifferentialEvo<double>>(bonds[i], irr_de);
		bonds[i].duration = bonds[i].macaulay_duration();
		std::cout << "Number of payements of bond " << i + 1 << " : " << bonds[i].cash_flows.size() << "\n";
		std::cout << "YTM of bond " << i + 1 << " : " << bonds[i].yield << "\n";
		std::cout << "Price of bond " << i + 1 << " : "  <<irr(bonds[i].yield, bonds[i].nominal_value, bonds[i].cash_flows, bonds[i].frequency) << "\n";
		std::cout << "Maculay Duration of bond " << i + 1 << " : " << bonds[i].duration << "\n";
	}
	std::cout << "Differential Evolution:" << '\n';
	DifferentialEvo<double> object_pricing(decision_variables, tol, 500, 0.5, 0.4, stdev);
	auto d = bond_pricing<double, DifferentialEvo<double>>(bonds, object_pricing);
	std::cout << "Optimum" << '\n';
	std::cout << d << " " << fitness_bond_pricing(d, bonds) << '\n';
	std::cout << "Particle Swarm:" << '\n';
	LocalBestPSO<double> object2(decision_variables, tol, 200, 1, 1, 2, 0.9, stdev);
	auto c = bond_pricing<double, LocalBestPSO<double>>(bonds, object2);
	std::cout << "Optimum" << '\n';
	std::cout << c << " " << fitness_bond_pricing(c, bonds) << '\n';
    return 0;
}

