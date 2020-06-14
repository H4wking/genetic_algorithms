#ifndef GENETIC_ALGORITHMS_CONFIG_HPP
#define GENETIC_ALGORITHMS_CONFIG_HPP

#define INT_BITS 10
#define DOUBLE_BITS 20
#define VARIABLES 2

#define POPULATION_SIZE 1000
#define GENERATIONS 300
#define MAX_NUM 512
#define MIN_NUM -512
#define MUTATION_CHANCE 0.1
#define THREADS 4

#define DOUBLE_MATE 1 // 1 - mate using split, 2 - mate using random genes method, 3 - mate using average method
#define BINARY_MATE 1 // 1 - mate using split, 2 - mate using random bits method

#define METHOD 2 // 1 - double, 2 - binary double,  3 - binary int

#define RUN_MPI false
#define MPI_FOR_DOUBLE true

#endif //GENETIC_ALGORITHMS_CONFIG_HPP
