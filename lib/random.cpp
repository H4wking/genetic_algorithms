#include <iostream>
#include "../inc/random.hpp"


/**
 * Generate random double from start number to end number
 *
 * @param start minimum number
 * @param end maximum number
 * @return random number
 */
double random_num(double start, double end) {
    double d = (double) rand() / RAND_MAX;
    return start + d * (end - start);
}


/**
 * Generate random integer from start number to end number
 *
 * @param start minimum number
 * @param end maximum number
 * @return random int
 */
int random_int(int start, int end) {
    int r = rand() % (end - start + 1);
    return start + r;
}


/**
 * Generate random binary number of length INT_BITS
 *
 * @return random binary string
 */
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


/**
 * Generate random binary number of length DOUBLE_BITS
 *
 * @return random binary string
 */
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


/**
 * Create vector of random binary numbers
 *
 * @return vector of string
 */
std::vector<std::string> create_int_ind() {
    std::vector<std::string> new_ind;
    for (int i = 0; i < VARIABLES; i++) {
        new_ind.push_back(random_bin_int());
    }
    return new_ind;
}


/**
 * Create vector of random binary numbers
 *
 * @return vector of string
 */
std::vector<std::string> create_double_ind() {
    std::vector<std::string> new_ind;
    for (int i = 0; i < VARIABLES; i++) {
        new_ind.push_back(random_bin_double());
    }
    return new_ind;
}


/**
 * Create vector of random doubles
 *
 * @return vector of doubles
 */
std::vector<double> create_idx_double() {
    std::vector<double> new_ind;
    for (int i = 0; i < VARIABLES; i++) {
        new_ind.push_back(mutated_gene());
    }
    return new_ind;
}


/**
 * Create random number from min to max
 *
 * @return double
 */
double mutated_gene() {
    return random_num(MIN_NUM, MAX_NUM);
}