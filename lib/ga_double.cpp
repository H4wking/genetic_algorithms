#include <iostream>
#include <vector>
#include <thread>
#include "../inc/ga_double.hpp"
#include "../inc/functions.hpp"
#include "../inc/random.hpp"
#include "../config.hpp"


/**
 * Create an Individual with parameters chromosome and func_res
 * chromosome is represented as a vector of doubles
 * Number of elements in vector depends on a number of dimension of functions that will be calculated
 * func_res - result of function calculated at given points - chromosome of Individual
 *
 */
Individual_double::Individual_double(std::vector<double> chromosome) {
    this->chromosome = chromosome;
    func_res = calculate_func_double();
}


/**
 * Create new individual by taking gene from both parents with equal chances and apply mutation
 *
 * @param par2 individual for mating
 * @return new Individual_double
 */
Individual_double Individual_double::mate_double(Individual_double par2) {
    std::vector<double> child_chromosome;

    double p = random_num(0, 1);
    int spl = random_int(0, VARIABLES - 2);

    if (p >= 0.5) {
        for (int i = 0; i <= spl; i++) {
            child_chromosome.push_back(chromosome[i]);
        }
        for (int j = spl + 1; j <= VARIABLES - 1; j++) {
            child_chromosome.push_back(par2.chromosome[j]);
        }
    } else {
        for (int i = 0; i <= spl; i++) {
            child_chromosome.push_back(par2.chromosome[i]);
        }
        for (int j = spl + 1; j <= VARIABLES - 1; j++) {
            child_chromosome.push_back(chromosome[j]);
        }
    }

    for (int i = 0; i < VARIABLES; i++) {
        p = random_num(0, 1);
        if (p <= MUTATION_CHANCE) {
            child_chromosome[i] = mutated_gene();
        }
    }

    return Individual_double(child_chromosome);
}



Individual_double Individual_double::mate2_double(Individual_double par2) {
    std::vector<double> child_chromosome;

    for (int i = 0; i < VARIABLES; i++) {
        double p = random_num(0, 1);
        if (p >= 0.5) {
            child_chromosome.push_back(chromosome[i]);
        } else {
            child_chromosome.push_back(par2.chromosome[i]);
        }
    }

    for (int i = 0; i < VARIABLES; i++) {
        double p = random_num(0, 1);
        if (p <= MUTATION_CHANCE) {
            child_chromosome[i] = mutated_gene();
        }
    }

    return Individual_double(child_chromosome);
}


/**
 * Create new individual by taking gene as average from both parents genes with and apply mutation
 *
 * @param par2 individual for mating
 * @return new Individual_double
 */
Individual_double Individual_double::mate_average(Individual_double par2) {
    std::vector<double> child_chromosome;

    for (int i = 0; i < VARIABLES; i++) {
        child_chromosome.push_back((chromosome[i] + par2.chromosome[i]) / 2);
    }

    for (int i = 0; i < VARIABLES; i++) {
        double p = random_num(0, 1);
        if (p <= MUTATION_CHANCE) {
            child_chromosome[i] = mutated_gene();
        }
    }

    return Individual_double(child_chromosome);
}


/**
 * Compare two Individuals by result of function.
 *
 * @relatesalso Individual_double
 */
bool operator<(const Individual_double &ind1, const Individual_double &ind2) {
    return ind1.func_res < ind2.func_res;
}


/**
 * Create initial population represented by a vector of Individuals
 *
 * @param pop_size number of individuals in initial generation
 * @return vector of Individuals
 */
std::vector<Individual_double> create_population_double(int pop_size) {
    std::vector<Individual_double> population;

    for(int i = 0; i < pop_size; i++)
    {
        std::vector<double> ind = create_idx_double();
        population.emplace_back(ind);
    }

    return population;
}


/**
 * Function for generating n offsprings from previous generation
 * Function is run in threads
 *
 * @param n number of offsprings
 * @param pop_size number of individuals in generation
 * @param prev vector of Individuals from previous generation
 * @return vector of Individuals of new generation
 */
void create_offsprings_thr(int n, int pop_size, std::vector<Individual_double> &prev_gen, std::vector<Individual_double> &new_gen) {
    for (int i = 0; i < n; i++) {
        int r = random_int(0, pop_size / 2);
        Individual_double parent1 = prev_gen[r];

        r = random_int(0, pop_size / 2);
        Individual_double parent2 = prev_gen[r];

        Individual_double offspring = parent1.mate2_double(parent2);

        new_gen.push_back(offspring);
    }
}


/**
 * Create new generation of Individuals
 * Add to new generation 10% of best individuals of previous population
 * and 90% of new individuals created from 50% of best from previous generation
 *
 * @param pop_size number of individuals in generation
 * @param prev vector of Individuals from previous generation
 * @return vector of Individuals of new generation
 */
std::vector<Individual_double> new_gen_double(int pop_size, std::vector<Individual_double> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual_double> new_generation;

    int s = (10 * pop_size) / 100;
    for(int i = 0; i < s; i++) {
        new_generation.push_back(prev[i]);
    }

    s = (90 * pop_size) / 100;

    std::vector<Individual_double> new_gens[THREADS];
    std::vector<std::thread> v;

    for (int i = 0; i < THREADS - 1; ++i)
        v.emplace_back(create_offsprings_thr, s / THREADS, pop_size,
                       std::ref(prev), std::ref(new_gens[i]));

    v.emplace_back(create_offsprings_thr, s - (s / THREADS) * (THREADS - 1), pop_size,
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


/**
 * Main function that runs genetic algorithm in next order:
 * - first create initial population of Individuals
 * - then create new generation of Individuals in loop given amount of times
 *
 *
 * @param pop_size number of individuals in generation
 * @param gen_num number of generations that will be created
 */
void run_ga_double(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_double> population = create_population_double(pop_size);

    while(generation <= gen_num)
    {
        population = new_gen_double(pop_size, population);

        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : population[0].chromosome) {
            std::cout << var << "; ";
        }
        std::cout << "\t";
        std::cout<< "Function result: "<< population[0].func_res << std::endl;

        generation++;
    }
}


/**
 * Calculate function for points that are represented as chromosome of Individual
 *
 * @return result of function
 */
double Individual_double::calculate_func_double() {
    return eggholder_function(chromosome[0], chromosome[1]);
}
