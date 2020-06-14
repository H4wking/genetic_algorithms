#include <iostream>
#include "../inc/ga_double_bit.hpp"


/**
 * Create an Individual with parameters chromosome and func_res
 * chromosome is represented as a vector of strings
 * Number of elements in vector depends on a number of dimension of functions that will be calculated
 * func_res - result of function calculated at given points - chromosome of Individual
 *
 */
Individual_bit::Individual_bit(std::vector<std::string> chromosome) {
    this->chromosome = chromosome;
    func_res = calculate_func_double_bit();
}


/**
 * Create new Individual as a result of mating of two Individuals.
 * Select random split point and create new individual by taking parts before and after split point from different parents.
 * Then apply mutation.
 *
 * @param par2 individual for mating
 * @return new individual
 */
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


/**
 * Create new individual by taking bits from both parents with equal chances and apply mutation
 *
 * @param par2 individual for mating
 * @return new individual
 */
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


/**
 * Compare two Individuals by result of function.
 *
 * @relatesalso Individual_bit
 */
bool operator<(const Individual_bit &ind1, const Individual_bit &ind2) {
    return ind1.func_res < ind2.func_res;
}


/**
 * Create initial population represented by a vector of Individuals
 *
 * @param pop_size number of individuals in initial generation
 * @return vector of Individuals
 */
std::vector<Individual_bit> create_population_double_bit(int pop_size) {
    std::vector<Individual_bit> population;

    for (int i = 0; i < pop_size; i++) {
        std::vector<std::string> ind = create_double_ind();
        population.emplace_back(ind);
    }

    return population;
}


/**
 * Function for generating n offsprings from previous generation
 * Function is run in threads
 * Does not return values, adds them to given vector
 *
 * @param n number of offsprings
 * @param pop_size number of individuals in generation
 * @param prev vector of Individuals from previous generation
 * @param vector for new generation
 */
void create_offsprings_thr_double_bit(int n, int pop_size, std::vector<Individual_bit> &prev_gen,
                                      std::vector<Individual_bit> &new_gen) {
    for (int i = 0; i < n; i++) {
        int r = random_int(0, pop_size / 2);
        Individual_bit parent1 = prev_gen[r];

        r = random_int(0, pop_size / 2);
        Individual_bit parent2 = prev_gen[r];

        if (BINARY_MATE == 1) {
            Individual_bit offspring = parent1.mate_double_bit(parent2);
            new_gen.push_back(offspring);
        } else {
            Individual_bit offspring = parent1.mate2_double_bit(parent2);
            new_gen.push_back(offspring);
        }
    }
}


/**
 * Create new generation of Individuals
 * Add to new generation 10% of best individuals of previous population
 * and 90% of new individuals created from 50% of best from previous generation
 * Creates threads that calculate new individuals and joins them afterwards
 *
 * @param pop_size number of individuals in generation
 * @param prev vector of Individuals from previous generation
 * @return vector of Individuals of new generation
 */
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


/**
 * Main function that runs genetic algorithm in next order:
 * - first create initial population of Individuals
 * - then create new generation of Individuals in loop given amount of times
 *
 *
 * @param pop_size number of individuals in generation
 * @param gen_num number of generations that will be created
 * @return vector of found minimum values
 */
std::vector<double> run_ga_double_bit(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_bit> population = create_population_double_bit(pop_size);

    std::vector<double> res;

    while (generation <= gen_num) {
        population = new_gen(pop_size, population);

        for (auto var : population[0].chromosome) {
            res.push_back(MIN_NUM + std::stoi(var, nullptr, 2) * ((MAX_NUM - MIN_NUM) / (pow(2, DOUBLE_BITS) - 1)));
        }

        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : res) {
            std::cout << var << "; ";
        }
        std::cout << "\t";
        std::cout << "Function result: " << population[0].func_res << std::endl;

        generation++;

        res.clear();
    }

    return res;
}


/**
 * Calculate function for points that are represented as chromosome of Individual
 *
 * @return result of function
 */
double Individual_bit::calculate_func_double_bit() {
    std::vector<double> chromosome_double;
    for (auto gene : chromosome) {
        double gene_double = MIN_NUM + std::stoi(gene, nullptr, 2) * ((MAX_NUM - MIN_NUM) / (pow(2, DOUBLE_BITS) - 1));
        chromosome_double.push_back(gene_double);
    }

    return eggholder_function(chromosome_double[0], chromosome_double[1]);
}


/**
 * Function for generating n offsprings from previous generation
 * Function is run in different processes using MPI
 *
 * @param n number of offsprings
 * @param prev vector of Individuals from previous generation
 */
std::vector<Individual_bit> create_offsprings_mpi_double_bit(int n, std::vector<Individual_bit> &prev_gen) {
    std::vector<Individual_bit> new_gen;

    for (int i = 0; i < n; i++) {
        int r = random_int(0, POPULATION_SIZE / 2 - 1);
        Individual_bit parent1 = prev_gen[r];

        r = random_int(0, POPULATION_SIZE / 2 - 1);
        Individual_bit parent2 = prev_gen[r];

        if (BINARY_MATE == 1) {
            Individual_bit offspring = parent1.mate_double_bit(parent2);
            new_gen.push_back(offspring);
        } else {
            Individual_bit offspring = parent1.mate2_double_bit(parent2);
            new_gen.push_back(offspring);
        }
    }

    return new_gen;
}


/**
 * Create new generation of Individuals
 * Add to new generation 10% of best individuals of previous population
 * and 90% of new individuals created from 50% of best from previous generation
 * Sends 50% of previous generation to other processes, which calculate and send back results
 * using mpi_new_gen_double() function
 *
 * @param pop_size number of individuals in generation
 * @param prev vector of Individuals from previous generation
 * @return vector of Individuals of new generation
 */
std::vector<Individual_bit> new_gen_mpi(int pop_size, std::vector<Individual_bit> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual_bit> new_generation;

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

    std::vector<Individual_bit> new_gens[numprocesses];

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
            MPI_Send(&str, DOUBLE_BITS + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        }
    }

    for (int i = 1; i <= numprocesses; i++) {
        std::vector<std::string> new_gen_str;

        std::string str;
        for (int j = 0; j < inds_process[i-1] * VARIABLES; j++) {
            MPI_Recv(&str, DOUBLE_BITS + 1, MPI_CHAR, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            new_gen_str.push_back(str);
        }

        for (int j = 0; j < inds_process[i-1] * VARIABLES; j = j + 2) {
            new_gens[i-1].push_back(Individual_bit(std::vector<std::string>{new_gen_str[j], new_gen_str[j+1]}));
        }
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
 * Creates new generation using MPI
 *
 * @param pop_size number of individuals in generation
 * @param gen_num number of generations that will be created
 * @return vector of found minimum values
 */
std::vector<double> run_ga_double_bit_mpi(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual_bit> population = create_population_double_bit(pop_size);

    std::vector<double> res;

    while (generation <= gen_num) {
        population = new_gen_mpi(pop_size, population);

        for (auto var : population[0].chromosome) {
            res.push_back(MIN_NUM + std::stoi(var, nullptr, 2) * ((MAX_NUM - MIN_NUM) / (pow(2, DOUBLE_BITS) - 1)));
        }

        std::cout << "Generation: " << generation << "\t";
        std::cout << "Minimum: ";
        for (auto var : res) {
            std::cout << var << "; ";
        }
        std::cout << "\t";
        std::cout << "Function result: " << population[0].func_res << std::endl;

        generation++;

        res.clear();
    }

    return res;
}


/**
 * Function for secondary processes to receive previous generation,
 * calculate new individuals and send them back
 */
void mpi_new_gen_double() {
    std::string str;
    int inds;

    std::vector<std::string> prev_str;
    std::vector<Individual_bit> prev;
    std::vector<Individual_bit> new_gen;
    std::vector<std::string> new_gen_str;

    for (int g = 0; g < GENERATIONS + 1; g++) {
        MPI_Recv(&inds, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < (POPULATION_SIZE / 2) * VARIABLES; i++) {
            MPI_Recv(&str, DOUBLE_BITS + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            prev_str.push_back(str);

        }

        for (int i = 0; i < (POPULATION_SIZE / 2) * VARIABLES; i = i + 2) {
            prev.push_back(Individual_bit(std::vector<std::string>{prev_str[i], prev_str[i + 1]}));
        }

        new_gen = create_offsprings_mpi_double_bit(inds, prev);

        for (int j = 0; j < inds; j++) {
            for (int i = 0; i < VARIABLES; i++) {
                new_gen_str.push_back(new_gen[j].chromosome[i]);
            }
        }

        for (auto s : new_gen_str) {
            MPI_Send(&s, DOUBLE_BITS + 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD);
        }

        prev_str.clear();
        prev.clear();
        new_gen.clear();
        new_gen_str.clear();
    }
}
