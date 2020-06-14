#include <iostream>
#include "../inc/functions.hpp"

//Main function for testing all methods
double eggholder_function(double x, double y) {
//    minimum: f(512,404.2319)=-959.6407
//    -512 < x,y < 512
    return -(y + 47) * sin(sqrt(abs(x / 2 + (y + 47))))
           - x * sin(sqrt(abs(x - (y + 47))));
}


double ackley_function(double x, double y) {
//    minimum: f(0,0)=0
//    -5 <= x,y <= 5
    return (-20 * exp(-0.2 * pow(0.5 * (pow(x, 2) + pow(y, 2)), 0.5))
            - exp(0.5 * (cos(2 * M_PI * x) + cos(2 * M_PI * y))) + exp(1) + 20);
}


double booth_function(double x, double y) {
    //    f(x,y)=(x+2y-7)^2+(2x+y-5)^{2} : f(1,3)=0
    //    -10 <= x,y <= 10
    return pow(x + 2 * y - 7, 2) + pow(2 * x + y - 5, 2);
}


double matyas_function(double x, double y) {
//    f(x,y)=0.26(x^2+y^2)-0.48xy
//    minium: f(0,0)=0
//    -10 <= x,y <= 10
    return 0.26 * (pow(x, 2) + pow(y, 2)) - 0.48 * x * y;
}


double rosenbrock_function (double x, double y) {
//   f(1.0,1.0)=0
//   -1.5 <= x, y <= 1.5
    return pow(1 - x, 2) + 100 * pow(y - pow(x, 2), 2);
}


double f2d_1(double x) {
//    2.7 <= x <= 7.5
    return sin(x) + sin(3.33 * x);
}


double f2d_2(double x) {
//    minimum f (0.9669) = -1.48907
//    0 <= x <= 1.2
    return -(1.4 - 3 * x) * sin(18 * x);
}
