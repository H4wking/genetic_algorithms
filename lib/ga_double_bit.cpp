#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <mutex>


#include "../inc/ga_double_bit.hpp"

#include "../inc/random.hpp"
#include "../config.hpp"



Individual_bit::Individual_bit(std::vector<std::string> chromosome) {
    this->chromosome = chromosome;
    func_res = calculate_func_double_bit();

}

Individual_bit Individual_bit::mate_double_bit(Individual_bit par2) {
    std::vector<std::string> child_chromosome;

    double p;
    int spl;
    std::string gene;

    for (int i = 0; i < VARIABLES; i++) {
        p = random_num(0, 1);
        spl = random_int(1, DOUBLE_BITS - 1);

        gene = "";
        if (p >= 0.5) {
            gene += chromosome[i].substr(0, spl);
            gene += par2.chromosome[i].substr(spl, DOUBLE_BITS);
        } else {
            gene += par2.chromosome[i].substr(0, spl);
            gene += chromosome[i].substr(spl, DOUBLE_BITS);
        }


        std::bitset<DOUBLE_BITS> gene_bin(gene);
        for (int j = 0; j < DOUBLE_BITS; j++) {
            p = random_num(0, 1);
            if (p <= MUTATION_CHANCE) {
                gene_bin.flip(j);
            }
        }

        child_chromosome.push_back(gene_bin.to_string());
    }


    return Individual_bit(child_chromosome);
}


Individual_bit Individual_bit::mate2_double_bit(Individual_bit par2) {
    std::vector<std::string> child_chromosome;

    double p;
    std::string gene;

    for (int i = 0; i < VARIABLES; i++) {
        gene = "";

        for (int j = 0; j < DOUBLE_BITS; j++) {
            p = random_num(0, 1);
            if (p >= 0.5) {
                gene += chromosome[i][j];
            } else {
                gene += par2.chromosome[i][j];
            }
        }

        std::bitset<DOUBLE_BITS> gene_bin(gene);
        for (int j = 0; j < DOUBLE_BITS; j++) {
            p = random_num(0, 1);
            if (p <= MUTATION_CHANCE) {
                gene_bin.flip(j);
            }
        }

        child_chromosome.push_back(gene_bin.to_string());
    }

    return Individual_bit(child_chromosome);
}

bool operator<(const Individual_bit &ind1, const Individual_bit &ind2) {
    return ind1.func_res < ind2.func_res;
}

std::vector<Individual_bit> create_population_double_bit(int pop_size) {
    std::vector<Individual_bit> population;

    for (int i = 0; i < pop_size; i++) {
        std::vector<std::string> ind = create_double_ind();
        population.emplace_back(ind);
    }

    return population;
}

void create_offsprings_thr_double_bit(int n, int pop_size, std::vector<Individual_bit> &prev_gen,
                                      std::vector<Individual_bit> &new_gen) {
    for (int i = 0; i < n; i++) {
        int r = random_int(0, pop_size / 2);
        Individual_bit parent1 = prev_gen[r];

        r = random_int(0, pop_size / 2);
        Individual_bit parent2 = prev_gen[r];

        Individual_bit offspring = parent1.mate_double_bit(parent2);

        new_gen.push_back(offspring);
    }
}

std::vector<Individual_bit> new_gen(int pop_size, std::vector<Individual_bit> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual_bit> new_generation;

    int s = (10 * pop_size) / 100;
    for (int i = 0; i < s; i++) {
        new_generation.push_back(prev[i]);
    }

    s = (90 * pop_size) / 100;

    std::vector<Individual_bit> new_gens[THREADS];
    std::vector<std::thread> v;

    for (int i = 0; i < THREADS - 1; ++i)
        v.emplace_back(create_offsprings_thr_double_bit, s / THREADS, pop_size,
                       std::ref(prev), std::ref(new_gens[i]));

    v.emplace_back(create_offsprings_thr_double_bit, s - (s / THREADS) * (THREADS - 1), pop_size,
                   std::ref(prev), std::ref(new_gens[THREADS - 1]));

    for (auto &t: v) {
        t.join();
    }

    new_generation.reserve(pop_size);
    for (auto g : new_gens) {
        new_generation.insert(new_generation.end(), g.begin(), g.end());
    }

    return new_generation;
}

void run_ga_double_bit(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_bit> population = create_population_double_bit(pop_size);

    while (generation <= gen_num) {
        population = new_gen(pop_size, population);


        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : population[0].chromosome) {
            std::cout << MIN_NUM + std::stoi(var, nullptr, 2) * ((MAX_NUM - MIN_NUM) / (pow(2, DOUBLE_BITS) - 1)) << "; ";
        }
        std::cout << "\t";
        std::cout << "Function result: " << population[0].func_res << std::endl;

        generation++;
    }
}


double Individual_bit::calculate_func_double_bit() {
    std::vector<double> chromosome_double;
    for (auto gene : chromosome) {
        double gene_double = MIN_NUM + std::stoi(gene, nullptr, 2) * ((MAX_NUM - MIN_NUM) / (pow(2, DOUBLE_BITS) - 1));
        chromosome_double.push_back(gene_double);
    }
    return -(chromosome_double[1] + 47) * sin(sqrt(abs(chromosome_double[0] / 2 + (chromosome_double[1] + 47))))
           - chromosome_double[0] * sin(sqrt(abs(chromosome_double[0] - (chromosome_double[1] + 47))));
}