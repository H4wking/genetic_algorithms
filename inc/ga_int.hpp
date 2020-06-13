#ifndef AKS_PROJECT_GA_INT_HPP
#define AKS_PROJECT_GA_INT_HPP


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

void run_ga_int(int gen_num, int pop_size);



#endif //AKS_PROJECT_GA_INT_HPP
