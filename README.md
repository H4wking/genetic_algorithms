# Parallel Library of Genetic Algorithms
## Project for Architecture of Computer Systems
__Task:__ create parallel library of genetic algorithms for function optimization  

Our library has three methods of genetic algorithms for function optimization:
* Genetic algorithm for double using bitwise crossover 
* Genetic algorithm for integer using bitwise crossover
* Genetic algorithm for double with three different crossover functions

First two approaches have two different crossover functions:
* Select random split point and create new individual by taking parts before and after split point from different parents. Then apply mutation.
* Create new individual by taking bits from both parents with equal chances and apply mutation.

Crossover functions of third approach:
* Select random split point and create new individual by taking parts before and after split point from different parents.
* Create new individual by taking gene from both parents with equal chances and apply mutation.
* Create new individual by taking gene as average from both parents genes and apply mutation.

Our library has two methods of parallelization:
* CPU using threads
* CPU using MPI

Each approach can be called using respective _run_ function.
For threads parallelization:
* run_ga_double_bit
* run_ga_int
* run_ga_double

For MPI parallelization:
* run_ga_double_bit_mpi
* run_ga_int_mpi

All of these functions have two parameters: number of generations and population size.

All parameters for our library are set in `config.hpp` file:
* INT_BITS - number of bits for integer method, sets interval
* DOUBLE_BITS - number of bits for double bit method, sets precision
* VARIABLES - number of variables in function
* POPULATION_SIZE - number of individuals in population
* GENERATIONS - number of generations
* MAX_NUM - upper boundary for double and double bit methods
* MIN_NUM - lower boundary for double and double bit methods
* MUTATION_CHANCE - probability of mutation
* THREADS - number of threads for thread parallelization
* METHOD - determine methos of genetic algorithm using thread parallelization: 1 - double, 2 - binary double,  3 - binary int
* RUN_MPI - boolean to choose thread or MPI parallelization
* MPI_FOR_DOUBLE - boolean to choose double bit or integer methods



## Description of functions using double bit method as example
Class Individual_bit representing individual which has two variables: chromosome - vector of strings containing binary numbers representing function variables, func_res - result of calculated function with variables from chromosome.

Methods: 
* 


