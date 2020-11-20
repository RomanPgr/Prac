#include <iostream>
#include "Windows.h"

double* prog(unsigned int N,
    const double* A,
    const double* C,
    const double* B,
    const double* F,
    const double &xi1,
    const double &xi2,
    const double &mu1,
    const double &mu2){

    double* alpha =  (double*) malloc(sizeof(double) * (N));
    double* beta =   (double*) malloc(sizeof(double) * (N));
    double* result = (double*) malloc(sizeof(double) * (N + 1));

    if (alpha == nullptr || beta == nullptr || result == nullptr){
        free(alpha);
        free(beta);
        free(result);
        return nullptr;
    }

    alpha[0] = xi1;
    beta[0] = mu1;

    for (unsigned int i = 0; i < N - 1; i++){
        alpha[i + 1] = B[i] / (C[i] - A[i] * alpha[i]);
        beta[i + 1] = (F[i] + A[i] * beta[i]) /
            (C[i] - A[i] * alpha[i]);
    }

    result[N] = (mu2 + xi2 * beta[N - 1]) / (1 - xi2 * alpha[N - 1]);

    for (int i = N; i > 0; i--){
        result[i - 1] = alpha[i - 1] * result[i] + beta[i - 1];
    }

    free(alpha);
    free(beta);

    return result;
}

int main()
{
    srand(174);

    for (unsigned int N = 100000000; N <= 200000000; N += 10000000){
        double* A = (double*) malloc(sizeof(double) * (N - 1));
        double* B = (double*) malloc(sizeof(double) * (N - 1));
        double* C = (double*) malloc(sizeof(double) * (N - 1));
        double* F = (double*) malloc(sizeof(double) * (N - 1));
        double* ideally = (double*) malloc(sizeof(double) * (N + 1));
        if (A == nullptr || B == nullptr || C == nullptr ||
            F == nullptr || ideally == nullptr){
            free(A);
            free(B);
            free(C);
            free(F);
            free(ideally);
            return -1;
        }
        double xi1, xi2, mu1, mu2;

        for (unsigned int i = 0; i < N + 1; i++){
            ideally[i] = (double) rand() / RAND_MAX;
        }

        for (unsigned int i = 0; i < N - 1; i++){
            A[i] = (double) rand() / RAND_MAX;
            B[i] = (double) rand() / RAND_MAX;
            C[i] = A[i] + B[i] + 1;
            F[i] = -A[i] * ideally[i] + C[i] * ideally[i + 1] - B[i] * ideally[i + 2];

        }

        mu1 = ideally[0];
        mu2 = ideally[N];
        xi1 = 0;
        xi2 = 0;

        double* result = prog(N, A, C, B, F, xi1, xi2, mu1, mu2);

        if (result == nullptr){
            free(A);
            free(B);
            free(C);
            free(F);
            free(ideally);
            return -1;
        }

        double max = 0, t;

        for (unsigned int i = 0; i < N + 1; i++){
            t = ideally[i] - result[i];
            if (t < 0){
                t = -t;
            }
            if (t > max){
                max = t;
            }
        }

        std::cout << "N = " << N << "\tmax = " << max << "\n";

        free(A);
        free(B);
        free(C);
        free(F);
        free(ideally);
        free(result);

        Sleep(10000);
    }

    std::cout << "\nSuccess" << std::endl;
    std::cin.get();

    return 0;
}
