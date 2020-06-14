#ifndef GENETIC_ALGORITHMS_RANDOM_HPP
#define GENETIC_ALGORITHMS_RANDOM_HPP

#include <vector>
#include <cstdlib>
#include <algorithm>
#include "../config.hpp"

double random_num(double start, double end);

int random_int(int start, int end);

std::string random_bin_int();

std::string random_bin_double();

std::vector<std::string> create_int_ind();

std::vector<std::string> create_double_ind();

double mutated_gene();


#endif //GENETIC_ALGORITHMS_RANDOM_HPP