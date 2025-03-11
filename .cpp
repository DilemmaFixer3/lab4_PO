#include <iostream>
#include <cmath>
#include <mpi.h>
#include <chrono>

double f(double x) {
    return sqrt(x - 1) + (1 / (x - 3));
}

int main(int argc, char* argv[]) {
    int rank, size;
    const int N = 1;
    const int a = N;
    const int b = N * 2;
    const double h = 0.1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int num_processes_to_run[] = { 1, 2, 4 };

    for (int i = 0; i < 3; ++i) {
        if (size == num_processes_to_run[i]) {
            double local_sum = 0.0;
            auto start_time = std::chrono::high_resolution_clock::now();

            for (int j = rank + 1; j <= b; j += size) {
                double x = a + (j - 1) * h;
                if (x <= b) {
                    local_sum += f(x);
                }
            }

            double total_sum;
            MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);

            if (rank == 0) {
                std::cout << "\nParallel computation with " << size << " processes:" << std::endl;
                std::cout << "Total sum: " << total_sum << std::endl;
                std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
            }
        }
    }

    MPI_Finalize();

    return 0;
}
