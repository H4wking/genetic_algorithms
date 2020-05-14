#include <iostream>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "../inc/random.hpp"

#include "../config.hpp"



double random_num(double start, double end) {
    double d = (double) rand() / RAND_MAX;
    return start + d * (end - start);
}

int random_int(int start, int end) {
    int r = rand() % (end - start + 1);
    return start + r;
}

std::string random_bin_int() {
    std::string bin = std::bitset<INT_BITS>(0).to_string();
    for (int i = 0; i < INT_BITS; i++) {
        double r = random_num(0, 1);
        if (r >= 0.5) {
            bin[i] = '1';
        }
    }

    return bin;
}

std::string random_bin_double() {
    std::string bin = std::bitset<DOUBLE_BITS>(0).to_string();
    for (int i = 0; i < DOUBLE_BITS; i++) {
        double r = random_num(0, 1);
        if (r >= 0.5) {
            bin[i] = '1';
        }
    }

    return bin;
}

std::vector<std::string> create_int_ind() {
    std::vector<std::string> new_ind;
    for (int i = 0; i < VARIABLES; i++) {
        new_ind.push_back(random_bin_int());
    }
    return new_ind;
}

std::vector<std::string> create_double_ind() {
    std::vector<std::string> new_ind;
    for (int i = 0; i < VARIABLES; i++) {
        new_ind.push_back(random_bin_double());
    }
    return new_ind;
}


double mutated_gene() {
    return random_num(MIN_NUM, MAX_NUM);
}