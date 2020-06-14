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




