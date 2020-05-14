#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "../inc/time.hpp"
#include <algorithm>
#include "../inc/ga_int.hpp"
#include "../inc/random.hpp"

#include "../config.hpp"

Individual_int::Individual_int(std::vector<std::string> chromosome) {
    this->chromosome = chromosome;
    func_res = calculate_func_int();
}

Individual_int Individual_int::mate_int(Individual_int par2) {
    std::vector<std::string> child_chromosome;
    std::string gene;
    int spl;
    double p;


    for (int i = 0; i < VARIABLES; i++) {
        p = random_num(0, 1);
        spl = random_int(0, VARIABLES - 1);

        gene = "";
        if (p >= 0.5) {
            gene += chromosome[i].substr(0, spl);
            gene += par2.chromosome[i].substr(spl, INT_BITS);

        } else {
            gene += par2.chromosome[i].substr(0, spl);
            gene += chromosome[i].substr(spl, INT_BITS);

        }
        std::bitset<INT_BITS> gene_bin(gene);
        for (int j = 0; j < INT_BITS; j++) {
            p = random_num(0, 1);
            if (p <= MUTATION_CHANCE) {
                gene_bin = gene_bin.flip(j);
            }
        }

        child_chromosome.push_back(gene_bin.to_string());
    }

    return Individual_int(child_chromosome);

}


bool operator<(const Individual_int &ind1, const Individual_int &ind2) {
    return ind1.func_res < ind2.func_res;
}

std::vector<Individual_int> create_population_int(int pop_size) {
    std::vector<Individual_int> population;

    for (int i = 0; i < pop_size; i++) {
        std::vector<std::string> ind = create_int_ind();

        population.emplace_back(ind);
    }

    return population;
}


std::vector<Individual_int> new_gen_int(int pop_size, std::vector<Individual_int> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual_int> new_generation;


    int s = (10 * pop_size) / 100;
    for (int i = 0; i < s; i++) {
        new_generation.push_back(prev[i]);
    }


    s = (90 * pop_size) / 100;


    for (int i = 0; i < s; i++) {

        int r = random_num(0, pop_size / 2);

        Individual_int parent1 = prev[r];
        r = random_num(0, pop_size / 2);
        Individual_int parent2 = prev[r];
        Individual_int offspring = parent1.mate_int(parent2);
        new_generation.push_back(offspring);

    }

    return new_generation;
}


void run_ga_int(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_int> population = create_population_int(pop_size);


    while (generation <= gen_num) {
        population = new_gen_int(pop_size, population);

        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : population[0].chromosome) {
            std::cout << MIN_NUM + std::stoi(var, nullptr, 2)  + 1<< "; ";
        }
        std::cout << "\t";
        std::cout << "Function result: " << population[0].func_res << std::endl;

        generation++;
    }


}


double Individual_int::calculate_func_int() {
    std::vector<int> chromosome_int;
    for (auto gene : chromosome) {
        chromosome_int.push_back(MIN_NUM + std::stoi(gene, nullptr, 2) + 1);
    }

    return -(chromosome_int[1] + 47) * sin(sqrt(abs(chromosome_int[0] / 2 + (chromosome_int[1] + 47))))
           - chromosome_int[0] * sin(sqrt(abs(chromosome_int[0] -  (chromosome_int[1] + 47))));

}
