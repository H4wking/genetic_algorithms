#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <cstdlib>

#define POPULATION_SIZE 5000
#define GENERATIONS 1000
#define MAX_NUM 512
#define MIN_NUM -512
#define MUTATION_CHANCE 0.5

double random_num(int start, int end) {
    double f = (double)rand() / RAND_MAX;
    return start + f * (end - start);
}

double mutated_gene() {
    return random_num(MIN_NUM, MAX_NUM);
}

std::pair<double, double> create()
{
    return std::pair<double, double>(mutated_gene(), mutated_gene());
}

class Individual {
public:
    std::pair<double, double> chromosome;
    double func_res;
    Individual(std::pair<double, double> chromosome);
    Individual mate(Individual parent2);
    double calculate_func();
};

Individual::Individual(std::pair<double, double> chromosome) {
    this->chromosome = chromosome;
    func_res = calculate_func();
}

Individual Individual::mate(Individual par2) {
    std::pair<double, double> child_chromosome;

    double p = random_num(0, 1);

    if (p >= 0.5) {
        child_chromosome.first = chromosome.first;
        child_chromosome.second = par2.chromosome.second;
    } else {
        child_chromosome.first = par2.chromosome.first;
        child_chromosome.second = chromosome.second;
    }


    p = random_num(0, 1);

    if (p <= MUTATION_CHANCE) {
        if (p <= MUTATION_CHANCE / 2) {
            child_chromosome.first = mutated_gene();
        } else {
            child_chromosome.second = mutated_gene();
        }
    }

    return Individual(child_chromosome);
}

bool operator<(const Individual &ind1, const Individual &ind2) {
    return ind1.func_res < ind2.func_res;
}

std::vector<Individual> create_population(int pop_size) {
    std::vector<Individual> population;

    for(int i = 0; i < pop_size; i++)
    {
        std::pair<double, double> ind = create();
        population.emplace_back(ind);
    }

    return population;
}

std::vector<Individual> new_gen(int pop_size, std::vector<Individual> prev) {
    sort(prev.begin(), prev.end());

    std::vector<Individual> new_generation;

    int s = (10 * pop_size) / 100;
    for(int i = 0; i < s; i++) {
        new_generation.push_back(prev[i]);
    }

    s = (90 * pop_size) / 100;
    for(int i = 0; i < s; i++)
    {
        int r = random_num(0, pop_size / 2);
        Individual parent1 = prev[r];
        r = random_num(0, pop_size / 2);
        Individual parent2 = prev[r];
        Individual offspring = parent1.mate(parent2);
        new_generation.push_back(offspring);
    }

    return new_generation;
}

void run(int gen_num, int pop_size) {
    int generation = 0;

    std::vector<Individual> population = create_population(pop_size);

    while(generation <= gen_num)
    {
        population = new_gen(pop_size, population);

        std::cout<< "Generation: " << generation << "\t";
        std::cout<< "(x, y): "<< population[0].chromosome.first << " " << population[0].chromosome.second <<"\t";
        std::cout<< "Function result: "<< population[0].func_res << std::endl;

        generation++;
    }
}

double Individual::calculate_func() {
    return -(chromosome.second + 47) * sin(sqrt(abs(chromosome.first / 2 + (chromosome.second + 47))))
           - chromosome.first * sin(sqrt(abs(chromosome.first -  (chromosome.second + 47))));
}

int main() {
    srand(time(nullptr));

    run(GENERATIONS, POPULATION_SIZE);
}
