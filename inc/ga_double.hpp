#ifndef GENETIC_ALGORITHMS_GA_DOUBLE_HPP
#define GENETIC_ALGORITHMS_GA_DOUBLE_HPP


std::vector<double> create_idx_double();

class Individual_double {
public:
    std::vector<double> chromosome;
    double func_res;
    Individual_double(std::vector<double> chromosome);
    Individual_double mate_double(Individual_double parent2);
    Individual_double mate2_double(Individual_double parent2);
    Individual_double mate_average(Individual_double parent2);
    double calculate_func_double();
};

bool operator<(const Individual_double &ind1, const Individual_double &ind2);

std::vector<Individual_double> create_population_double(int pop_size);

void create_offsprings_thr(int n, int pop_size, std::vector<Individual_double> &prev_gen, std::vector<Individual_double> &new_gen);

std::vector<Individual_double> new_gen_double(int pop_size, std::vector<Individual_double> prev);

void run_ga_double(int gen_num, int pop_size);



#endif //GENETIC_ALGORITHMS_GA_DOUBLE_HPP
