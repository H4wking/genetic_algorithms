#ifndef GENETIC_ALGORITHMS_GA_DOUBLE_BIT_HPP
#define GENETIC_ALGORITHMS_GA_DOUBLE_BIT_HPP

#include <vector>
#include <cmath>
#include <thread>
#include <mutex>
#include <mpi.h>
#include "../inc/functions.hpp"
#include "../inc/random.hpp"
#include "../config.hpp"

class Individual_bit {
public:
    std::vector<std::string> chromosome;
    double func_res;

    Individual_bit(std::vector<std::string> chromosome);

    Individual_bit mate_double_bit(Individual_bit parent2);

    Individual_bit mate2_double_bit(Individual_bit parent2);

    double calculate_func_double_bit();
};

bool operator<(const Individual_bit &ind1, const Individual_bit &ind2);

std::vector<Individual_bit> create_population_double_bit(int pop_size);

void create_offsprings_thr_double_bit(int n, int pop_size, std::vector<Individual_bit> &prev_gen,
                                      std::vector<Individual_bit> &new_gen);

std::vector<Individual_bit> new_gen(int pop_size, std::vector<Individual_bit> prev);

std::vector<double> run_ga_double_bit(int gen_num, int pop_size);

std::vector<Individual_bit> create_offsprings_mpi_double_bit(int n, std::vector<Individual_bit> &prev_gen);

std::vector<Individual_bit> new_gen_mpi(int pop_size, std::vector<Individual_bit> prev);

std::vector<double> run_ga_double_bit_mpi(int gen_num, int pop_size);

void mpi_new_gen_double();

#endif //GENETIC_ALGORITHMS_GA_DOUBLE_BIT_HPP