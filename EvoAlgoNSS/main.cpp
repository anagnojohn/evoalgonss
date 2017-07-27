// ConsoleApplication1.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "svensson.h"
#include "bond.h"
#include "irr.h"
#include "bondpricing.h"
#include "dependencies.h"
#include "benchmarks.h"
#include <array>

int main()
{
	Bond<double> bond1 { 0.06, 101.657, 100, 2, "2016-03-30", "2019-06-24" };
	Bond<double> bond2 { 0.09, 99.675, 100, 2, "2016-03-30", "2022-11-22" };
	Bond<double> bond3 { 0.07, 126.888, 100, 2, "2016-03-30", "2024-01-15" };
	Bond<double> bond4 { 0.04, 136.2, 100, 2, "2016-03-30", "2027-07-09" };
	Bond<double> bond5 { 0.08, 100.22, 100, 2, "2016-03-30", "2030-04-28" };
	Bond<double> bond6 { 0.05, 125.98, 100, 2, "2016-03-30", "2039-03-15" };
	Bond<double> bond7 { 0.06, 132.88, 100, 2, "2016-03-30", "2043-03-15" };
	Bond<double> bond8 {0.11, 150.33, 100, 2, "2016-03-30", "2048-07-28" };
	std::vector<Bond<double>> bonds { bond1, bond2, bond3, bond4, bond5, bond6, bond7, bond8};
	GAstruct<double> ga_irr { 0.4, 0.35, 200, 0.0001, 200 };
	GAstruct<double> ga_pricing{ 0.4, 0.35, 200, 0.0000001, 200 };
	//PSOstruct<double> pso_pricing { 1, 1, 5, { 1000, 1000, 1000, 1000 }, 2, 200, 0.001, 200 };
	DEstruct<double> de_pricing { 0.6, 1, 200, 0.0001, 200 };
	DEstruct<double> de_irr { 0.6, 1, 200, 0.0000001, 200 };
	benchmarkcurvefitting(bonds, de_irr, de_pricing);
	//benchmarkcurvefitting(bonds, ga_irr, ga_pricing);
    return 0;
}