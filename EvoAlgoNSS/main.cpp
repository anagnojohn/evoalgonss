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
		decision_variables[i][0]= b0;
		decision_variables[i][1] = b1;
		decision_variables[i][2] = b2;
		decision_variables[i][3] = b3;
		decision_variables[i][4] = tau1;
		decision_variables[i][5] = tau2;
	}
	//std::cout << bond1.settlement_date;
	std::vector<double> bond_yields = { 0.043355724, 0.055688708, 0.040845871 ,0.02564443,
		0.057006077, 0.037201461 ,0.041785886 ,0.049908247 };
	std::vector<double> duration = { 2.944988754, 4.966178711, 6.279674883, 9.474865358,
		8.416273259, 14.93089635, 15.38446779, 12.22184684 };
	auto f = [&](const std::vector<double>& solution) { return fitness_svensson<double>(solution, bond_yields, duration); };
	std::vector<double> stdev(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		stdev[j] = 0.7;
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
	
	Bond<double> bond1(0.06, 1.0166, 1, 2, "2016-03-30", "2019-06-24");
	Bond<double> bond2(0.09, 0.9968, 1, 2, "2016-03-30", "2022-11-22");
	Bond<double> bond3(0.07, 1.2689, 1, 2, "2016-03-30", "2024-01-15");
	Bond<double> bond4(0.04, 1.3620, 1, 2, "2016-03-30", "2027-07-09");
	Bond<double> bond5(0.08, 1.0022, 1, 2, "2016-03-30", "2030-04-28");
	Bond<double> bond6(0.05, 1.2598, 1, 2, "2016-03-30", "2039-03-15");
	Bond<double> bond7(0.06, 1.3288, 1, 2, "2016-03-30", "2043-03-15");
	Bond<double> bond8(0.11, 1.5033, 1, 2, "2016-03-30", "2048-07-28");
	bonds.push_back(bond1);
	bonds.push_back(bond2);
	bonds.push_back(bond3);
	bonds.push_back(bond4);
	bonds.push_back(bond5);
	bonds.push_back(bond6);
	bonds.push_back(bond7);
	bonds.push_back(bond8);
	Bond<double> bondtest(0.06, 1001, 1000, 2, "2016-03-30", "2019-06-24");
	std::cout << bondtest.duration;
	auto f_irr = [&](const std::vector<double>& solution) { return fitness_irr(solution, bondtest.price, bondtest.nominal_value, bondtest.cash_flows, bondtest.frequency);};
	size_t ndv_irr = 1;
	std::vector<double> stdev_irr(1);
	for (auto j = 0; j < ndv_irr; ++j)
	{
		stdev_irr[j] = 0.7;
	}
	std::cout << "Differential Evolution:" << '\n';
	auto f_pricing = [&](const std::vector<double>& solution) { return fitness_bond_pricing<double>(solution, bonds); };
	DifferentialEvo<double, decltype(f_pricing)> object_pricing(decision_variables, tol, 0.0, 500, 0.5, 0.4, stdev);
	std::vector<double> d = object_pricing.differential_evo(f_pricing);
	std::cout << "Optimum" << '\n';
	std::cout << d << " " << f_pricing(d) << '\n';
	std::cout << "Particle Swarm:" << '\n';
	LocalBestPSO<double, decltype(f_pricing)> object2(decision_variables, tol, 0.0, 200, 1, 1, 2, 0.9, stdev);
	auto c = object2.lbest_pso(f_pricing);
	std::cout << "Optimum" << '\n';
	std::cout << c << " " << f(c) << '\n';
	DifferentialEvo<double, decltype(f_irr)> irr_de(init_rate, tol, 0.0, 200, 0.5, 0.4, stdev_irr);
	//std::vector<double> res = irr_de.differential_evo(f_irr);
	//bond1.yield = res[0];
	//std::cout << bond1.yield;
    return 0;
}

