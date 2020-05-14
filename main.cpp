#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include <ctime>
#include <cstdlib>
#include <thread>
#include <mutex>
#include "inc/time.hpp"
#include <algorithm>


#include "inc/random.hpp"

#include "inc/ga_int.hpp"
#include "inc/ga_double.hpp"

#include "inc/ga_double_bit.hpp"



#include "config.hpp"

int main() {
    srand(time(nullptr));

    auto start = get_current_time_fenced();

//    run_ga_double(GENERATIONS, POPULATION_SIZE);

    run_ga_double_bit(GENERATIONS, POPULATION_SIZE);

//    run_ga_int(GENERATIONS, POPULATION_SIZE);

    auto finish = get_current_time_fenced();

    auto total_time  = finish - start;

    std::cout << to_us(total_time) / 1000 << std::endl;
}
