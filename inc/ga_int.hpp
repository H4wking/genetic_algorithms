#ifndef GENETIC_ALGORITHMS_GA_INT_HPP
#define GENETIC_ALGORITHMS_GA_INT_HPP

#include <vector>
#include <thread>
#include <algorithm>
#include <mpi.h>
#include "../inc/random.hpp"
#include "../config.hpp"
#include "../inc/functions.hpp"

class Individual_int {
public:
    std::vector<std::string> chromosome;
    double func_res;

    Individual_int(std::vector<std::string> chromosome);

    Individual_int mate_int(Individual_int parent2);

    Individual_int mate2_int(Individual_int parent2);

    double calculate_func_int();
};

bool operator<(const Individual_int &ind1, const Individual_int &ind2);

std::vector<Individual_int> create_population_int(int pop_size);

void create_offsprings_thr_int(int n, int pop_size, std::vector<Individual_int> &prev_gen,
                               std::vector<Individual_int> &new_gen);

std::vector<Individual_int> new_gen_int(int pop_size, std::vector<Individual_int> prev) ;

std::vector<int> run_ga_int(int gen_num, int pop_size);

std::vector<Individual_int> create_offsprings_mpi_int(int n, std::vector<Individual_int> &prev_gen);

std::vector<Individual_int> new_gen_int_mpi(int pop_size, std::vector<Individual_int> prev) ;

std::vector<int> run_ga_int_mpi(int gen_num, int pop_size);

void mpi_new_gen_int();


#endif //GENETIC_ALGORITHMS_GA_INT_HPP