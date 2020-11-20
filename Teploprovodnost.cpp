#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>

double* prog(unsigned int N,
    const double* A,
    const double* C,
    const double* B,
    const double* F,
    const double xi1,
    const double xi2,
    const double mu1,
    const double mu2){

    double* alpha = (double*)malloc(sizeof(double) * (N));
    double* beta = (double*)malloc(sizeof(double) * (N));
    double* result = (double*)malloc(sizeof(double) * (N + 1));

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


double* temp(unsigned int N,
    double(*f)(double x, double t),
    double(*k)(double x),
    double(*mu_1)(double t),
    double(*mu_2)(double t),
    double(*fi)(double x),
    double sigma,
    double h,
    double tau,
    double t_max){

    unsigned int T = (int)(t_max / tau) + 1;

    double* A = (double*)malloc(sizeof(double) * (N - 1));
    double* B = (double*)malloc(sizeof(double) * (N - 1));
    double* C = (double*)malloc(sizeof(double) * (N - 1));
    double* F = (double*)malloc(sizeof(double) * (N - 1));
    double* result = (double*)malloc(sizeof(double) * (T * (N + 1)));

    if (result == nullptr || A == nullptr || B == nullptr || C == nullptr || F == nullptr){
        free(A);
        free(B);
        free(C);
        free(F);
        free(result);
        return nullptr;
    }

    {
        double tmp = sigma * tau / (h * h);
        result[0] = (*mu_1)(0);
        result[N] = (*mu_2)(0);
        for (unsigned int i = 0; i < N - 1; i++){
            A[i] = tmp * (*k)(h * (i + 0.5));
            B[i] = tmp * (*k)(h * (i + 1.5));
            C[i] = A[i] + B[i] + 1.0;
            result[i + 1] = (*fi)((i + 1) * h);
        }
    }
    
    for (unsigned int j = 1; j < T; j++){
        for (unsigned int i = 0; i < N - 1; i++){
            F[i] =  tau * (sigma * (*f)((i + 1) * h, (j + 0) * tau) + (1 - sigma) * (*f)((i + 1) * h, (j - 1) * tau)) +
                tau * (1 - sigma) / (h * h) * ((
                (*k)(i + 1.5) * (result[(N + 1) * (j - 1) + (i + 2)] - result[(N + 1) * (j - 1) + i + 1]) - 
                (*k)(i + 0.5) * (result[(N + 1) * (j - 1) + (i + 1)] - result[(N + 1) * (j - 1) + i + 0]))) +
                result[(N + 1) * (j - 1) + i + 1];
        }
        double* tempr = prog(N, A, C, B, F, 0.0, 0.0, (*mu_1)(j * tau), (*mu_2)(j * tau));
        if (tempr == nullptr){
            free(A);
            free(B);
            free(C);
            free(F);
            free(result);
            return nullptr;
        }
        for (unsigned int i = 0; i < N + 1; i++){
            result[j * (N + 1) + i] = tempr[i];
        }

        free(tempr);
    }

    free(A);
    free(B);
    free(C);
    free(F);

    return result;
}

double k(double x){
    return 1.0;
}

double f(double x, double t){
    return 0.0;
}

double mu_1(double t){
    return 0.0;
}

double mu_2(double t){
    return 0.0;
}

double fi(double x){
    double l = 1;
    return std::sin(M_PI * l * x);
}

double resh(double x, double t){
    double l = 1;
    return (std::pow(M_E, -M_PI * M_PI * l * l * t) * std::sin(M_PI * l * x));
}

int main()
{
    int N = 100;
    double sigma = 1.0,
        h = 0.01,
        tau =  0.01,
        T_max = 3,
        max = 0,
        gr;
    int T = (int)((T_max / tau) + 1);

    double* yt = temp(N, f, k, mu_1, mu_2, fi, sigma, h, tau, T_max);
    if (yt == nullptr){
        return -1;
    }
    std::cout << "[";
    unsigned int vt = T;// / 100;
    for (int i = 0; i < vt; i++){
        std::cout /*<< i << " " << N + 1 */<< "[";
        for (int j = 0; j < N + 1; j++){
            std::cout << yt[i * (N + 1) + j];
            if (j != N){
                std::cout << ", ";
            }
            gr = yt[i * (N + 1) + j] - resh(j * h, i * tau);
            if (gr < 0){
                gr = -gr;
            }
            if (gr > max){
                max = gr;
            }
        }
        std::cout << "]";
        if (i != vt - 1){
            std::cout << ", ";
        }
    }
    std::cout << "]";
    std::cout << "\nmax = " << max << "\n";
    //std::cin.get();
    free(yt);

    return 0;
}
