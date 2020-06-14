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
#include <mpi.h>



#include "inc/random.hpp"

#include "inc/ga_int.hpp"
#include "inc/ga_double.hpp"
#include "inc/ga_double_bit.hpp"


#include "config.hpp"

int main(int argc, char *argv[]) {
    if (!RUN_MPI) {
        srand(time(nullptr));

        auto start = get_current_time_fenced();

//        run_ga_double(GENERATIONS, POPULATION_SIZE);

        run_ga_double_bit(GENERATIONS, POPULATION_SIZE);

//        run_ga_int(GENERATIONS, POPULATION_SIZE);

        auto finish = get_current_time_fenced();

        auto total_time = finish - start;

        std::cout << to_us(total_time) / 1000 << std::endl;
    } else {
        int rank;

        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        if (rank == 0) {
            srand(time(nullptr));

            auto start = get_current_time_fenced();

            if (MPI_FOR_DOUBLE) {
                run_ga_double_bit_mpi(GENERATIONS, POPULATION_SIZE);
            } else {
                run_ga_int_mpi(GENERATIONS, POPULATION_SIZE);
            }

            auto finish = get_current_time_fenced();

            auto total_time = finish - start;

            std::cout << to_us(total_time) / 1000 << std::endl;
        } else {
            if (MPI_FOR_DOUBLE) {
                mpi_new_gen_double();
            } else {
                mpi_new_gen_int();
            }
        }

        MPI_Finalize();
    }
}
