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
    const double h = 0.000001;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Послідовне обчислення
    double total_sum_sequential = 0.0;
    auto start_time_sequential = std::chrono::high_resolution_clock::now();

    for (int i = 1; i <= b; ++i) {
        double x = a + (i - 1) * h;
        if (x <= b) {
            total_sum_sequential += f(x);
        }
    }

    auto end_time_sequential = std::chrono::high_resolution_clock::now();
    auto duration_sequential = std::chrono::duration_cast<std::chrono::microseconds>(end_time_sequential - start_time_sequential);

    if (rank == 0) {
        std::cout << "Sequential computation:" << std::endl;
        std::cout << "Total sum: " << total_sum_sequential << std::endl;
        std::cout << "Time taken: " << duration_sequential.count() << " microseconds" << std::endl;
    }

    // Паралельне обчислення
    double local_sum_parallel = 0.0;
    auto start_time_parallel = std::chrono::high_resolution_clock::now();

    for (int i = rank + 1; i <= b; i += size) {
        double x = a + (i - 1) * h;
        if (x <= b) {
            local_sum_parallel += f(x);
        }
    }

    double total_sum_parallel;
    MPI_Reduce(&local_sum_parallel, &total_sum_parallel, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    auto end_time_parallel = std::chrono::high_resolution_clock::now();
    auto duration_parallel = std::chrono::duration_cast<std::chrono::microseconds>(end_time_parallel - start_time_parallel);

    if (rank == 0) {
        std::cout << "\nParallel computation:" << std::endl;
        std::cout << "Total sum: " << total_sum_parallel << std::endl;
        std::cout << "Time taken: " << duration_parallel.count() << " microseconds" << std::endl;
    }

    MPI_Finalize();

    return 0;
}


//#include <iostream>
//#include <cmath>
//#include <mpi.h>
//#include <chrono>
//
//// Функція, яка обчислює значення функції f(x)
//double f(double x) {
//    return sqrt(x - 1) + (1 / (x - 3));
//}
//
//int main(int argc, char* argv[]) {
//    int rank, size;
//    const int N = 1; // Змінна для встановлення меж інтегрування
//    const int a = N; // Нижня межа інтегрування
//    const int b = N * 2;  // Верхня межа інтегрування
//    const double h = 0.000001; // Крок інтегрування
//
//    // Ініціалізація MPI
//    MPI_Init(&argc, &argv);
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Отримання номера процесу
//    MPI_Comm_size(MPI_COMM_WORLD, &size); // Отримання загальної кількості процесів
//
//    int num_processes_to_run[] = { 1, 2, 4 }; // Масив, що містить кількість процесів, які будуть виконувати обчислення
//
//    // Паралельне обчислення для різної кількості процесів
//    for (int i = 0; i < 3; ++i) {
//        if (size == num_processes_to_run[i]) {
//            double local_sum = 0.0; // Локальна сума, яка буде обчислена кожним процесом
//            auto start_time = std::chrono::high_resolution_clock::now(); //початок відліку часу
//
//            // Цикл для обчислення частини інтегралу, яка припадає на кожний процес
//            for (int j = rank + 1; j <= b; j += size) {
//                double x = a + (j - 1) * h; // Визначення значення x
//                if (x <= b) {
//                    local_sum += f(x); // Додавання до локальної суми значення функції f(x)
//                }
//            }
//
//            double total_sum;
//            MPI_Reduce(&local_sum, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);// Зведення локальних сум у загальну суму
//
//            auto end_time = std::chrono::high_resolution_clock::now();// Кінцева мітка часу
//            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time); // Обчислення тривалості виконання
//
//            if (rank == 0) {
//                // Виведення результатів для головного процесу
//                std::cout << "\nParallel computation with " << size << " processes:" << std::endl;
//                std::cout << "Total sum: " << total_sum << std::endl;
//                std::cout << "Time taken: " << duration.count() << " microseconds" << std::endl;
//            }
//        }
//    }
//
//    MPI_Finalize(); // Завершення роботи з MPI
//
//    return 0;
//}



//#include <iostream> 
//#include <mpi.h> 
//using namespace std;
//
//int main(int argc, char* argv[]) {
//	int num_procs, rank; 
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &num_procs); 
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	const int n = 70; 
//	const int m = 120; 
//	const int p = 200;
//
//	int A[n][m], B[m][p], C[n][p], Cb[n][p];
//
//	// Initialize matrices A and B on rank 0 
//	if (rank == 0) {
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < m; j++) {
//				A[i][j] = i + j;
//			}
//		}
//
//		for (int i = 0; i < m; i++) {
//			for (int j = 0; j < p; j++) {
//				B[i][j] = i + j;
//			}
//		}
//
//		// Output matrix A
//		cout << "Matrix A:" << endl; 
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < m; j++) {
//				cout << A[i][j] << " ";
//			}
//			cout << endl;
//		}
//		cout << endl;
//
//		// Output matrix B
//		cout << "Matrix B:" << endl; 
//		for (int i = 0; i < m; i++) {
//			for (int j = 0; j < p; j++) {
//				cout << B[i][j] << " ";
//			}
//			cout << endl;
//		}
//		cout << endl;
//	}
//
//	// Broadcast matrices A and B to all processes
//	MPI_Bcast(A, n* m, MPI_INT, 0, MPI_COMM_WORLD); 
//	MPI_Bcast(B, m* p, MPI_INT, 0, MPI_COMM_WORLD);
//
//	// Calculate rows and columns per process 
//	int rows_per_process = n / num_procs; 
//	int cols_per_process = p / num_procs;
//
//	// Matrix multiplication
//	for (int i = rank * rows_per_process; i < (rank + 1) * rows_per_process; i++) {
//		for (int j = 0; j < p; j++) {
//			C[i][j] = 0;
//			for (int k = 0; k < m; k++) {
//				C[i][j] += A[i][k] * B[k][j];
//			}
//
//		}
//	}
//
//	// Gather results
//	MPI_Reduce(C, Cb, n * p, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//
//	// Output results from rank 0 
//	if (rank == 0) {
//		cout << "Matrix Cb:" << endl; 
//		for (int i = 0; i < n; i++) {
//			for (int j = 0; j < p; j++) {
//				cout << Cb[i][j] << " ";
//			}
//			cout << endl;
//		}
//	}
//
//	MPI_Finalize(); 
//	return 0;
//}

//#include <math.h>
//#include "mpi.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <iostream>
//
//using namespace std;
//
//int ProcRank = 0;
//int ProcNum = 0;
//double* pMatrixX;
//double* pMatrixy1;
//double* pMatrixy2;
//double* pMatrixY3;
//int Size = 10;
//int Size2 = 40;
//
//void RandomDataInitialization(double* pMatrix, double* pVector, int Size) {
//    int i, j;
//    for (i = 0; i < Size; i++) {
//        pVector[i] = rand() / double(100);
//        for (j = 0; j < Size2; j++)
//            pMatrix[i * Size2 + j] = rand() / double(100);
//    }
//}
//
//void ResultReplication(double* pProcResult, double* pResult, int Size, int RowNum) {
//    int* pReceiveNum;
//    int* pReceiveInd;
//    int RestRows = Size2;
//    int i;
//    pReceiveNum = new int[ProcNum];
//    pReceiveInd = new int[ProcNum];
//    pReceiveInd[0] = 0;
//    pReceiveNum[0] = Size2 / ProcNum;
//    for (i = 1; i < ProcNum; i++) {
//        RestRows -= pReceiveNum[i - 1];
//        pReceiveNum[i] = RestRows / (ProcNum - i);
//        pReceiveInd[i] = pReceiveInd[i - 1] + pReceiveNum[i - 1];
//    }
//    MPI_Allgatherv(pProcResult, pReceiveNum[ProcRank], MPI_DOUBLE, pResult, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
//    delete[] pReceiveNum;
//    delete[] pReceiveInd;
//}
//
//void ParallelResultCalculation(double* pProcRows, double* pVector, double* pProcResult, int Size, int RowNum) {
//    int i, j;
//    for (i = 0; i < RowNum; i++) {
//        pProcResult[i] = 0;
//        for (j = 0; j < Size2; j++)
//            pProcResult[i] += pProcRows[i * Size2 + j] * pVector[j];
//    }
//}
//
//void ParallelResultCalculationM(double* pProcRows, double* pProcRows2, double* pProcResult, int Size, int RowNum) {
//    int i, j;
//    for (i = 0; i < RowNum; i++) {
//        double summ = 0;
//        pProcResult[i] = 0;
//        for (j = 0; j < Size2; j++)
//            summ += pProcRows[i * Size2 + j] * pProcRows2[j];
//        pProcResult[i] = summ;
//    }
//}
//
//void ParallelResultCalculationS(double* pProcRows, double* pProcRows2, double* pProcResult, int Size, int RowNum) {
//    int i, j;
//    for (i = 0; i < RowNum; i++) {
//        double summ = 0;
//        pProcResult[i] = 0;
//        for (j = 0; j < Size2; j++)
//            summ += pProcRows[i * Size2 + j] + pProcRows2[j];
//        pProcResult[i] = summ;
//    }
//}
//
//void ProcessTermination(double* pMatrix, double* pVector, double* pResult, double* pProcRows, double* pProcResult) {
//    if (ProcRank == 0)
//        delete[] pMatrix;
//    delete[] pVector;
//    delete[] pResult;
//    delete[] pProcRows;
//    delete[] pProcResult;
//}
//
//void DataDistribution(double* pMatrix, double* pProcRows, double* pVector, int Size, int ColNum) {
//    int* pSendNum;
//    int* pSendInd;
//    int RestCol = Size2;
//    MPI_Bcast(pVector, Size2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    pSendInd = new int[ProcNum];
//    pSendNum = new int[ProcNum];
//    ColNum = (Size2 / ProcNum);
//    pSendNum[0] = ColNum * Size2;
//    pSendInd[0] = 0;
//    for (int i = 1; i < ProcNum; i++) {
//        RestCol -= ColNum;
//        ColNum = RestCol / (ProcNum - i);
//        pSendNum[i] = ColNum * Size2;
//        pSendInd[i] = pSendInd[i - 1] + pSendNum[i - 1];
//    }
//    MPI_Scatterv(pMatrix, pSendNum, pSendInd, MPI_DOUBLE, pProcRows, pSendNum[ProcRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    delete[] pSendNum;
//    delete[] pSendInd;
//}
//
//void PrintMatrix(double* pMatrix, int RowCount, int ColCount) {
//    int i, j;
//    for (i = 0; i < RowCount; i++) {
//        for (j = 0; j < ColCount; j++)
//            printf("%7.4f ", pMatrix[i * ColCount + j]);
//        printf("\n");
//    }
//}
//
//void PrintVector(double* pVector, int Size) {
//    int i;
//    for (i = 0; i < Size; i++)
//        printf("%7.4f\n", pVector[i]);
//}
//
//void ProcessInitialization(double*& pMatrix, double*& pVector, double*& pResult, double*& pProcRows, double*& pProcResult, int& Size, int& RowNum) {
//    int RestRows;
//    int i;
//    MPI_Bcast(&Size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    RestRows = Size;
//    for (i = 0; i < ProcRank; i++)
//        RestRows = RestRows - RestRows / (ProcNum - i);
//    RowNum = RestRows / (ProcNum - ProcRank);
//    pVector = new double[Size2];
//    pResult = new double[Size];
//    pProcRows = new double[RowNum * Size2];
//    pProcResult = new double[RowNum];
//    if (ProcRank == 0) {
//        pMatrix = new double[Size * Size2];
//        RandomDataInitialization(pMatrix, pVector, Size2);
//    }
//}
//
//int main(int argc, char* argv[]) {
//    double* pMatrixA;
//    double* pMatrixA1;
//    double* pMatrixA2;
//    double* pMatrixB2;
//    double* pVectorb1;
//    double* pVectorc1;
//    double* pVector;
//    double* pResult;
//    double* pProcRowsA;
//    double* pProcRowsA1;
//    double* pProcRowsA2;
//    double* pProcRowsB2;
//    double* pProcResult;
//    int RowNum;
//    double Start, Finish, Duration;
//    MPI_Init(&argc, &argv);
//    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
//    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
//    if (ProcRank == 0) {
//        printf("Size: %dx%d\nProcs: %d\n", Size, Size2, ProcNum);
//    }
//    double* pVectorb = new double[Size2];
//    ProcessInitialization(pMatrixA, pVector, pResult, pProcRowsA, pProcResult, Size2, RowNum);
//    ProcessInitialization(pMatrixA1, pVectorc1, pResult, pProcRowsA1, pProcResult, Size2, RowNum);
//    ProcessInitialization(pMatrixA2, pVectorb1, pResult, pProcRowsA2, pProcResult, Size2, RowNum);
//    ProcessInitialization(pMatrixB2, pVectorc1, pResult, pProcRowsB2, pProcResult, Size2, RowNum);
//    if (ProcRank == 0)
//        Start = MPI::Wtime();
//    DataDistribution(pMatrixA, pProcRowsA, pVector, Size2, RowNum);
//    DataDistribution(pMatrixA1, pProcRowsA1, pVector, Size2, RowNum);
//    DataDistribution(pMatrixA2, pProcRowsA2, pVector, Size2, RowNum);
//    //DataDistribution(pMatrixB2, pProcRows
//{
