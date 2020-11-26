#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>

double resh(double x, double t);

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


double temp(unsigned int N,
    unsigned int S_max,
    unsigned int chastota,
    double(*f)(double x, double t),
    double(*k)(double x),
    double(*mu_1)(double t),
    double(*mu_2)(double t),
    double(*fi)(double x),
    double h,
    double tau,
    double t_max){

    unsigned int T = (int)(t_max / tau);// +1;

    double* A = (double*)malloc(sizeof(double) * (N - 1));
    double* B = (double*)malloc(sizeof(double) * (N - 1));
    double* C = (double*)malloc(sizeof(double) * (N - 1));
    double* F = (double*)malloc(sizeof(double) * (N - 1));
    double* result_time = (double*)malloc(sizeof(double) * (N + 1));
    double* result_iter = (double*)malloc(sizeof(double) * (N + 1));

    double max = 0.0, modul;

    if (result_iter == nullptr || result_time == nullptr || A == nullptr || B == nullptr || C == nullptr || F == nullptr){
        free(A);
        free(B);
        free(C);
        free(F);
        free(result_iter);
        free(result_time);
    }
    
    for (unsigned int i = 0; i < N + 1; i++){
        result_iter[i] = (*fi)(h * i);
        std::cout << result_iter[i] << ' ';
    }
    std::cout << '\n';

    memcpy(result_time, result_iter, sizeof(double) * (N + 1));
    double tmp = tau / (h * h);

    for (unsigned int j = 1; j < T; j++){
        for (unsigned s = 1; s <= S_max; s++){
            for (unsigned int i = 0; i < N - 1; i++){
                A[i] = tmp * (*k)((result_iter[i + 0] + result_iter[i + 1]) / 2);
                B[i] = tmp * (*k)((result_iter[i + 1] + result_iter[i + 2]) / 2);
                C[i] = A[i] + B[i] + 1.0;
                F[i] = result_time[i + 1] + tau * ((*f)((i + 1) * h, j * tau));
            }

            double* tempr = prog(N, A, C, B, F, 0.0, 0.0, (*mu_1)(j * tau), (*mu_2)(j * tau));
            if (tempr == nullptr){
                free(A);
                free(B);
                free(C);
                free(F);
                free(result_iter);
                free(result_time);
            }

            memcpy(result_iter, tempr, sizeof(double) * (N + 1));
            free(tempr);
        }
        memcpy(result_time, result_iter, sizeof(double) * (N + 1));

        for (unsigned int i = 0; i < N + 1; i++){
            if ((modul = std::abs(result_time[i] - resh(i * h, j * tau))) > max){
                max = modul;
            }
        }

        if (j % chastota == 0){
            for (unsigned int v = 0; v < N + 1; v++){
                std::cout << result_time[v] << ' ';
            }
            std::cout << '\n';
        }
    }

    std::cout.flush();

    free(A);
    free(B);
    free(C);
    free(F);
    free(result_iter);
    free(result_time);
    
    return max;
}

double k(double u){
    double a = 1, 
        b = 3, 
        s = 2;
    return a + b * std::pow(u, s);
}

double f(double x, double t){
    //return std::pow(x, 10);
    return 1.0;
}

double mu_1(double t){
    return 0.0;
}

double mu_2(double t){
    return 0.0;
}

double fi(double x){
    double l = 1.0;
    return std::sin(M_PI * l * x);
}

double resh(double x, double t){
    double l = 1;
    return (std::pow(M_E, -M_PI * M_PI * l * l * t) * std::sin(M_PI * l * x));
}

int main()
{
    unsigned int N = 100, S_max = 1, chastota = 1;
    double h = 0.01,
        tau = 0.1,
        T_max = 3.0;
    int T = (int)((T_max / tau) + 1);

    double m = temp(N, S_max, chastota, f, k, mu_1, mu_2, fi, h, tau, T_max);

    //std::cout << "max = " << m << std::endl;
    //std::cin.get();

    return 0;
}
