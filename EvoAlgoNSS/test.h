#pragma once

double dc1_max = 5.0;
double dc1_min = -5.0;
double dc2_max = 5.0;
double dc2_min = -5.0;
double dc3_max = 5.0;
double dc3_min = -5.0;
size_t npop = 1200;
std::default_random_engine generator;
std::uniform_real_distribution<> distribution(0.0, 1.0);
std::vector<std::vector<double>> decision_variables;
for (auto i = 0; i < npop; ++i)
{
	double x = dc1_min + (distribution(generator)) * (dc1_max - dc1_min);
	double y = dc2_min + (distribution(generator)) * (dc2_max - dc2_min);
	double z = dc3_min + (distribution(generator)) * (dc3_max - dc3_min);
	decision_variables.push_back({ x,y,z });
	//decision_variables.push_back({ x,y });
}
auto f = [&](std::vector<double> dv) { return std::pow(dv[0], 2) + std::pow(dv[1], 2) + 1 + std::pow(dv[2], 2) + 1; };
//auto f = [&](std::vector<double> dv) { return std::pow(dv[0], 2) + std::pow(dv[1], 2) + 1; };
//std::cout << "Decision Variables: " << std::endl;
//for (const auto& p : decision_variables)
//{
//	std::cout << p << " " << f(p) << std::endl;
//}