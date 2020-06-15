#ifndef GENETIC_ALGORITHMS_FUNCTIONS_HPP
#define GENETIC_ALGORITHMS_FUNCTIONS_HPP

#include <cmath>
#include <vector>

double function_for_optimization(std::vector<double> variables);

double eggholder_function(std::vector<double> variables);

double ackley_function(std::vector<double> variables);

double booth_function(std::vector<double> variables);

double matyas_function(std::vector<double> variables);

double rosenbrock_function (std::vector<double> variables);

double f2d_1(std::vector<double> variables);

double f2d_2(std::vector<double> variables);


#endif