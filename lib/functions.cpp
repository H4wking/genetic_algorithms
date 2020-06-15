#include <iostream>
#include "../inc/functions.hpp"

// Main function for testing all methods, can have any formula
double function_for_optimization(std::vector<double> variables) {
    return -(variables[1] + 47) * sin(sqrt(abs(variables[0] / 2 + (variables[1] + 47))))
           - variables[0] * sin(sqrt(abs(variables[0] - (variables[1] + 47))));
}


double eggholder_function(std::vector<double> variables) {
//    minimum: f(512,404.2319)=-959.6407
//    -512 < x,y < 512
    return -(variables[1] + 47) * sin(sqrt(abs(variables[0] / 2 + (variables[1] + 47))))
           - variables[0] * sin(sqrt(abs(variables[0] - (variables[1] + 47))));
}


double ackley_function(std::vector<double> variables) {
//    minimum: f(0,0)=0
//    -5 <= x,y <= 5
    return (-20 * exp(-0.2 * pow(0.5 * (pow(variables[0], 2) + pow(variables[1], 2)), 0.5))
            - exp(0.5 * (cos(2 * M_PI * variables[0]) + cos(2 * M_PI * variables[1]))) + exp(1) + 20);
}


double booth_function(std::vector<double> variables) {
    //    f(x,y)=(x+2y-7)^2+(2x+y-5)^{2} : f(1,3)=0
    //    -10 <= x,y <= 10
    return pow(variables[0] + 2 * variables[1] - 7, 2) + pow(2 * variables[0] + variables[1] - 5, 2);
}


double matyas_function(std::vector<double> variables) {
//    f(x,y)=0.26(x^2+y^2)-0.48xy
//    minium: f(0,0)=0
//    -10 <= x,y <= 10
    return 0.26 * (pow(variables[0], 2) + pow(variables[1], 2)) - 0.48 * variables[0] * variables[1];
}


double rosenbrock_function (std::vector<double> variables) {
//   f(1.0,1.0)=0
//   -1.5 <= x, y <= 1.5
    return pow(1 - variables[0], 2) + 100 * pow(variables[1] - pow(variables[0], 2), 2);
}


double f2d_1(std::vector<double> variables) {
//    2.7 <= x <= 7.5
    return sin(variables[0]) + sin(3.33 * variables[0]);
}


double f2d_2(std::vector<double> variables) {
//    minimum f (0.9669) = -1.48907
//    0 <= x <= 1.2
    return -(1.4 - 3 * variables[0]) * sin(18 * variables[0]);
}
