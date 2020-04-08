#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>

#define POPULATION_SIZE 10000
#define GENERATIONS 1000
#define MAX_NUM 5
#define MIN_NUM -5
#define MUTATION_CHANCE 0.4

double random_num(int start, int end)
{
    double f = (double)rand() / RAND_MAX;
    return start + f * (end - start);
}

double mutated_genes()
{
    return random_num(MIN_NUM, MAX_NUM);
}

std::pair<double, double> create_gnome()
{
    return std::pair<double, double>(mutated_genes(), mutated_genes());
}

class Individual
{
public:
    std::pair<double, double> chromosome;
    double fitness;
    Individual(std::pair<double, double> chromosome);
    Individual mate(Individual parent2);
    double cal_fitness();
};

Individual::Individual(std::pair<double, double> chromosome)
{
    this->chromosome = chromosome;
    fitness = cal_fitness();
};

Individual Individual::mate(Individual par2)
{
    std::pair<double, double> child_chromosome;

    double p = random_num(0, 1);

    if (p >= 0.5) {
        child_chromosome.first = chromosome.first;
        child_chromosome.second = par2.chromosome.second;
    } else {
        child_chromosome.first = par2.chromosome.first;
        child_chromosome.second = chromosome.second;
    }


    double p2 = random_num(0, 1);

    if (p2 <= MUTATION_CHANCE) {
        if (p2 <= MUTATION_CHANCE / 2) {
            child_chromosome.first = mutated_genes();
        } else {
            child_chromosome.second = mutated_genes();
        }
    }

    return Individual(child_chromosome);
};

double Individual::cal_fitness()
{
    return pow(1 - chromosome.first, 2) + 100 * pow(chromosome.second - pow(chromosome.first, 2), 2);
};

bool operator<(const Individual &ind1, const Individual &ind2)
{
    return ind1.fitness < ind2.fitness;
}

int main()
{
    srand((unsigned)(time(0)));

    int generation = 0;

    std::vector<Individual> population;

    for(int i = 0; i < POPULATION_SIZE; i++)
    {
        std::pair<double, double> gnome = create_gnome();
        population.push_back(Individual(gnome));
    }

    while(generation <= GENERATIONS)
    {
        sort(population.begin(), population.end());

        std::vector<Individual> new_generation;

        int s = (10 * POPULATION_SIZE) / 100;
        for(int i = 0; i < s; i++) {
            new_generation.push_back(population[i]);
        }

        s = (90 * POPULATION_SIZE) / 100;
        for(int i = 0; i < s; i++)
        {
            int r = random_num(0, POPULATION_SIZE / 2);
            Individual parent1 = population[r];
            r = random_num(0, POPULATION_SIZE / 2);
            Individual parent2 = population[r];
            Individual offspring = parent1.mate(parent2);
            new_generation.push_back(offspring);
        }
        population = new_generation;
        std::cout<< "Generation: " << generation << "\t";
        std::cout<< "(x, y): "<< population[0].chromosome.first << " " << population[0].chromosome.second <<"\t";
        std::cout<< "Fitness: "<< population[0].fitness << std::endl;

        generation++;
    }
}
