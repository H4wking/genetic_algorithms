#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <cstdlib>

#define POPULATION_SIZE 1000
#define GENERATIONS 1000
#define MAX_NUM 512
#define MIN_NUM -512
#define MUTATION_CHANCE 0.7
#define VARIABLES 2

double random_num(int start, int end) {
    double d = (double)rand() / RAND_MAX;
    return start + d * (end - start);
}

int random_int(int start, int end) {
    int r = rand() % (end - start + 1);
    return start + r;
}

double mutated_gene() {
    return random_num(MIN_NUM, MAX_NUM);
}

std::vector<double> create() {
    std::vector<double> new_ind;
    for (int i = 0; i < VARIABLES; i++) {
        new_ind.push_back(mutated_gene());
    }
    return new_ind;
}

class Individual {
public:
    std::vector<double> chromosome;
    double func_res;
    Individual(std::vector<double> chromosome);
    Individual mate(Individual parent2);
    Individual mate2(Individual parent2);
    Individual mate_average(Individual parent2);
    double calculate_func();
};

Individual::Individual(std::vector<double> chromosome) {
    this->chromosome = chromosome;
    func_res = calculate_func();
}

Individual Individual::mate(Individual par2) {
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

    return Individual(child_chromosome);
}

Individual Individual::mate2(Individual par2) {
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

    return Individual(child_chromosome);
}

Individual Individual::mate_average(Individual par2) {
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

    return Individual(child_chromosome);
}

bool operator<(const Individual &ind1, const Individual &ind2) {
    return ind1.func_res < ind2.func_res;
}

std::vector<Individual> create_population(int pop_size) {
    std::vector<Individual> population;

    for(int i = 0; i < pop_size; i++)
    {
        std::vector<double> ind = create();
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
        int r = random_int(0, pop_size / 2);
        Individual parent1 = prev[r];
        r = random_int(0, pop_size / 2);
        Individual parent2 = prev[r];
        Individual offspring = parent1.mate2(parent2);
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

double Individual::calculate_func() {
    return -(chromosome[1] + 47) * sin(sqrt(abs(chromosome[0] / 2 + (chromosome[1] + 47))))
           - chromosome[0] * sin(sqrt(abs(chromosome[0] -  (chromosome[1] + 47))));
}

int main() {
    srand(time(nullptr));

    run(GENERATIONS, POPULATION_SIZE);
}
