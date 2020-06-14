#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <algorithm>
#include <mpi.h>
#include "../inc/time.hpp"
#include "../inc/ga_int.hpp"
#include "../inc/random.hpp"
#include "../config.hpp"
#include "../inc/functions.hpp"

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


Individual_int Individual_int::mate2_int(Individual_int par2) {
    std::vector<std::string> child_chromosome;

    double p;
    std::string gene;

    for (int i = 0; i < VARIABLES; i++) {
        gene = "";

        for (int j = 0; j < INT_BITS; j++) {
            p = random_num(0, 1);
            if (p >= 0.5) {
                gene += chromosome[i][j];
            } else {
                gene += par2.chromosome[i][j];
            }
        }

        std::bitset<INT_BITS> gene_bin(gene);
        for (int j = 0; j < INT_BITS; j++) {
            p = random_num(0, 1);
            if (p <= MUTATION_CHANCE) {
                gene_bin.flip(j);
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


void create_offsprings_thr_int(int n, int pop_size, std::vector<Individual_int> &prev_gen,
                               std::vector<Individual_int> &new_gen) {
    for (int i = 0; i < n; i++) {
        int r = random_int(0, pop_size / 2);
        Individual_int parent1 = prev_gen[r];

        r = random_int(0, pop_size / 2);
        Individual_int parent2 = prev_gen[r];

        Individual_int offspring = parent1.mate2_int(parent2);

        new_gen.push_back(offspring);
    }
}


std::vector<Individual_int> new_gen_int(int pop_size, std::vector<Individual_int> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual_int> new_generation;


    int s = (10 * pop_size) / 100;
    for (int i = 0; i < s; i++) {
        new_generation.push_back(prev[i]);
    }


    s = (90 * pop_size) / 100;

    std::vector<Individual_int> new_gens[THREADS];
    std::vector<std::thread> v;

    for (int i = 0; i < THREADS - 1; ++i)
        v.emplace_back(create_offsprings_thr_int, s / THREADS, pop_size,
                       std::ref(prev), std::ref(new_gens[i]));

    v.emplace_back(create_offsprings_thr_int, s - (s / THREADS) * (THREADS - 1), pop_size,
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


void run_ga_int(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_int> population = create_population_int(pop_size);


    while (generation <= gen_num) {
        population = new_gen_int(pop_size, population);

        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : population[0].chromosome) {
            std::cout << MIN_NUM + std::stoi(var, nullptr, 2) + 1 << "; ";
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

    return eggholder_function(chromosome_int[0], chromosome_int[1]);
}

std::vector<Individual_int> create_offsprings_mpi_int(int n, std::vector<Individual_int> &prev_gen) {
    std::vector<Individual_int> new_gen;

    for (int i = 0; i < n; i++) {
        int r = random_int(0, POPULATION_SIZE / 2 - 1);
        Individual_int parent1 = prev_gen[r];

        r = random_int(0, POPULATION_SIZE / 2 - 1);
        Individual_int parent2 = prev_gen[r];

        Individual_int offspring = parent1.mate2_int(parent2);

        new_gen.push_back(offspring);
    }

    return new_gen;
}


std::vector<Individual_int> new_gen_int_mpi(int pop_size, std::vector<Individual_int> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual_int> new_generation;

    int s = (10 * pop_size) / 100;
    for (int i = 0; i < s; i++) {
        new_generation.push_back(prev[i]);
    }

    s = (90 * pop_size) / 100;

    int numprocesses;

    MPI_Comm_size(MPI_COMM_WORLD, &numprocesses);

    numprocesses--;

    int avinds = s / numprocesses;
    int extra = s % numprocesses;

    int inds;

    std::vector<Individual_int> new_gens[numprocesses];

    std::vector<int> inds_process;

    std::vector<std::string> prev_str;

    for (int j = 0; j < pop_size / 2; j++) {
        for (int i = 0; i < VARIABLES; i++) {
            prev_str.push_back(prev[j].chromosome[i]);
        }
    }

    for (int i = 1; i <= numprocesses; i++) {
        inds = (i <= extra) ? avinds + 1 : avinds;
        inds_process.push_back(inds);

        MPI_Send(&inds, 1, MPI_INT, i, 0, MPI_COMM_WORLD);

        for (auto str : prev_str) {
            MPI_Send(&str, INT_BITS + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
    }

    for (int i = 1; i <= numprocesses; i++) {
        std::vector<std::string> new_gen_str;

        std::string str;
        for (int j = 0; j < inds_process[i-1] * VARIABLES; j++) {
            MPI_Recv(&str, INT_BITS + 1, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            new_gen_str.push_back(str);
        }

        for (int j = 0; j < inds_process[i-1] * VARIABLES; j = j + 2) {
            new_gens[i-1].push_back(Individual_int(std::vector<std::string>{new_gen_str[j], new_gen_str[j+1]}));
        }
    }

    new_generation.reserve(pop_size);
    for (auto g : new_gens) {
        new_generation.insert(new_generation.end(), g.begin(), g.end());
    }

    return new_generation;
}


void run_ga_int_mpi(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_int> population = create_population_int(pop_size);


    while (generation <= gen_num) {
        population = new_gen_int_mpi(pop_size, population);

        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : population[0].chromosome) {
            std::cout << MIN_NUM + std::stoi(var, nullptr, 2) + 1 << "; ";
        }
        std::cout << "\t";
        std::cout << "Function result: " << population[0].func_res << std::endl;

        generation++;
    }


}

void mpi_new_gen_int() {
    std::string str;
    int inds;

    std::vector<std::string> prev_str;
    std::vector<Individual_int> prev;
    std::vector<Individual_int> new_gen;
    std::vector<std::string> new_gen_str;

    for (int g = 0; g < GENERATIONS + 1; g++) {
        MPI_Recv(&inds, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < (POPULATION_SIZE / 2) * VARIABLES; i++) {
            MPI_Recv(&str, INT_BITS + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            prev_str.push_back(str);

        }

        for (int i = 0; i < (POPULATION_SIZE / 2) * VARIABLES; i = i + 2) {
            prev.push_back(Individual_int(std::vector<std::string>{prev_str[i], prev_str[i + 1]}));
        }

        new_gen = create_offsprings_mpi_int(inds, prev);

        for (int j = 0; j < inds; j++) {
            for (int i = 0; i < VARIABLES; i++) {
                new_gen_str.push_back(new_gen[j].chromosome[i]);
            }
        }

        for (auto s : new_gen_str) {
            MPI_Send(&s, INT_BITS + 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
        }

        prev_str.clear();
        prev.clear();
        new_gen.clear();
        new_gen_str.clear();
    }
}
