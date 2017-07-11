// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"
#include "svensson.h"
#include "datehandler.h"

int main()
{
	test();
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
		

	std::vector<double> bond_yields = { 0.043355724, 0.055688708, 0.040845871 ,0.02564443,
		0.057006077, 0.037201461 ,0.041785886 ,0.049908247 };
	std::vector<double> duration = { 2.944988754, 4.966178711, 6.279674883, 9.474865358,
		8.416273259, 14.93089635, 15.38446779, 12.22184684 };
	auto f = [&](const std::vector<double>& solution) { return fitness_svensson<double>(solution, bond_yields, svensson<double>, duration); };
	std::vector<double> stdev(ndv);
	for (auto j = 0; j < ndv; ++j)
	{
		stdev[j] = 0.7;
	}
	std::cout << "Genetic Algorithm:" << '\n';
	GeneticAlgo<double, decltype(f)> object1(decision_variables, tol, 0.0, 200, 0.4, 0.35, stdev);
	std::cout << "Optimum" << '\n';
	//auto b = object1.genetic_algo(f);
	//std::cout << b << " " << f(b) << '\n';
	std::cout << "Particle Swarm:" << '\n';
	LocalBestPSO<double, decltype(f)> object2(decision_variables, tol, 0.0, 200, 1, 1, 2, 0.9, stdev);
	auto c = object2.lbest_pso(f);
	std::cout << "Optimum" << '\n';
	std::cout << c << " " << f(c) << '\n';
	std::cout << "Differential Evolution:" << '\n';
	DifferentialEvo<double, decltype(f)> object3(decision_variables, tol, 0.0, 200, 0.5, 0.4, stdev);
	auto d = object3.differential_evo(f);
	std::cout << "Optimum" << '\n';
	std::cout << d << " " << f(d) << '\n';
    return 0;
}

