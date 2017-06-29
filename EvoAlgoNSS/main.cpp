// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "geneticalgo.h"
#include "local_best_pso.h"
#include "differentialevo.h"
#include "svensson.h"

int main()
{
	double b0 = 4.336 / 100;
	double b1 = 4.991 / 100 - b0;
	double b2 = 0;
	double b3 = 0;
	double b4 = 0;
	double tau1 = 1;
	double tau2 = 2;
	size_t npop = 500;
	size_t ndv = 6;
	double maturity = 3;
	double tol = 0.001;
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
	std::vector<double> bond_yields = { 0.004336, 0.005569, 0.004085, 0.002564, 0.005701, 0.003720, 0.004179, 0.004991 };
	auto f = [&](const std::vector<double>& solution) { return fitness_svensson<double>(solution, bond_yields, svensson<double>, maturity); };
	std::cout << "Genetic Algorithm:" << '\n';
	//auto b = genetic_algo<double>(decision_variables, f, tol, 0.0, 200);
	//std::cout << "Optimum" << '\n';
	//std::cout << b << " " << f(b) << '\n';
	std::cout << "Particle Swarm:" << '\n';
	auto c = lbest_pso<double>(decision_variables, f, tol, 0.0, 200);
	std::cout << "Optimum" << '\n';
	std::cout << c << " " << f(c) << '\n';
	std::cout << "Differential Evolution:" << '\n';
	auto d = differential_evo<double>(decision_variables, f, tol, 0.0, 200);
	std::cout << "Optimum" << '\n';
	std::cout << d << " " << f(d) << '\n';
    return 0;
}

