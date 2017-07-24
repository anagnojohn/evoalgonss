// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"
#include "svensson.h"
#include "bond.h"
#include "irr.h"
#include "bondpricing.h"
#include "dependencies.h"

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
	double maturity = 3;
	const double tol = 0.0001;
	std::vector<double> decision_variables = { b0, b1, b2, b3, tau1, tau2 };
	const size_t ndv = decision_variables.size();
	size_t npop_irr = 100;
	std::vector<double> init_rate = { 0.2 };
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
	EAparams<double> irr_params(init_rate, npop_irr, 0.000000001, 200, stdev_irr);
	DifferentialEvo<double> irr_de(0.5, 0.9);
	for (auto i = 0; i < bonds.size(); ++i)
	{
		bonds[i].yield = setyield<double, DifferentialEvo<double>>(bonds[i], irr_de, irr_params);
		bonds[i].duration = bonds[i].macaulay_duration();
		std::cout << "Number of payements of bond " << i + 1 << " : " << bonds[i].cash_flows.size() << "\n";
		std::cout << "YTM of bond " << i + 1 << " : " << bonds[i].yield << "\n";
		std::cout << "Price of bond " << i + 1 << " : "  <<irr(bonds[i].yield, bonds[i].nominal_value, bonds[i].cash_flows, bonds[i].frequency) << "\n";
		std::cout << "Maculay Duration of bond " << i + 1 << " : " << bonds[i].duration << "\n";
	}
	std::cout << "Differential Evolution:" << '\n';
	EAparams<double> pricing_params(decision_variables, npop, tol, 500, stdev);
	//DifferentialEvo<double> de_pricing(0.5, 0.4);
	//auto d = bond_pricing<double, DifferentialEvo<double>>(bonds, de_pricing, pricing_params);
	//std::cout << "Optimum" << '\n';
	//std::cout << d << " " << fitness_bond_pricing(d, bonds) << '\n';
	std::cout << "Genetic Algorithm:" << '\n';
	GeneticAlgo<double> object1(0.4, 0.35);
	auto b = bond_pricing<double, GeneticAlgo<double>>(bonds, object1, pricing_params);
	std::cout << "Optimum" << '\n';
	std::cout << b << " " << fitness_bond_pricing(b, bonds) << '\n';
	std::cout << "Particle Swarm:" << '\n';
	std::vector<double> vmax(ndv);
	for (auto& p : vmax)
	{
		p = 1000;
	}
	LocalBestPSO<double> object2(1, 1, 5, vmax, 2);
	auto c = bond_pricing<double, LocalBestPSO<double>>(bonds, object2, pricing_params);
	std::cout << "Optimum" << '\n';
	std::cout << c << " " << fitness_bond_pricing(c, bonds) << '\n';
    return 0;
}

