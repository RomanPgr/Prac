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


int koleb(
    unsigned int N,
    double *pred,
    double *cur,
    double(*f)(double x, double t),
    double(*mu_1)(double t),
    double(*mu_2)(double t),
    double sigma,
    double h,
    double tau,
    unsigned int T,
    unsigned int period){

    double max = 0.0;

    for (unsigned int j = 0; j <= N; j++){
        std::cout << pred[j] << ' ';
    }
    std::cout << '\n';

    double* A = (double*)malloc(sizeof(double) * (N - 1));
    double* B = (double*)malloc(sizeof(double) * (N - 1));
    double* C = (double*)malloc(sizeof(double) * (N - 1));
    double* F = (double*)malloc(sizeof(double) * (N - 1));

    if (A == nullptr || B == nullptr || C == nullptr || F == nullptr){
        free(A);
        free(B);
        free(C);
        free(F);
        return -1;
    }

    {
        double tmp = sigma * (tau * tau) / (h * h);
        for (unsigned int i = 0; i < N - 1; i++){
            A[i] = tmp;
            B[i] = tmp;
            C[i] = tmp + tmp + 1.0;
        }
    }

    for (unsigned int j = 2; j < T; j++){
        for (unsigned int i = 0; i < N - 1; i++){
            F[i] =
                2 * cur[i + 1] -
                pred[i + 1] +
                tau * tau * ((*f)((i + 1) * h, (j - 1) * tau) +
                (sigma / (h * h)) * (pred[i] - 2 * pred[i + 1] + pred[i + 2]) +
                ((1 - 2 * sigma) / (h * h)) * (cur[i] - 2 * cur[i + 1] + cur[i + 2]));
        }
        double* result = prog(N, A, C, B, F, 0.0, 0.0, (*mu_1)(j * tau), (*mu_2)(j * tau));
        if (result == nullptr){
            free(A);
            free(B);
            free(C);
            free(F);
            free(result);
            return -1;
        }

        for (unsigned int m = 0; m <= N; m++){
            if (std::abs(cur[m] - resh(m * h, (j - 1) * tau)) > max){
                max = std::abs(cur[m] - resh(m * h, (j - 1) * tau));
            }
        }


        free(pred);
        pred = cur;
        cur = result;

        if (j % ((int)(T / period)) == 0){
            for (unsigned int m = 0; m <= N; m++){
                std::cout << cur[m] << ' ';
            }
            std::cout << '\n';
        }
    }
    //std::cout << max;

    std::cout.flush();
    free(A);
    free(B);
    free(C);
    free(F);
    free(pred);
    free(cur);

    return 0;
}

double f(double x, double t){
    return 0.0;
    return (- x * x * std::sin(t));
}

double mu_1(double t){
    double x = 0.0;
    return std::cos(M_PI * t) * std::cos(M_PI * x);
    return 2 * t + (x * x - 2) * std::sin(t) + std::cos(M_PI * t) * std::cos(M_PI * x);
    return 0.0;
}

double mu_2(double t){
    double x = 1.0;
    return std::cos(M_PI * t) * std::cos(M_PI * x);
    return 2 * t + (x * x - 2) * std::sin(t) + std::cos(M_PI * t) * std::cos(M_PI * x);
    return 0.0;// std::sin(M_PI * 1.5);
}

double fi(double x){
    //double l = 1.5;
    //return std::sin(M_PI * l * x);
    //return x * x;
    return std::cos(M_PI * x);
}

double psi(double x){
    //double l = 1.5;
    //return std::sin(M_PI * l * x) * 0.9;
    return 0.0;
}

double resh(double x, double t){
    return std::cos(M_PI * t) * std::cos(M_PI * x);
    return (2 * t + (x * x - 2) * std::sin(t) + std::cos(M_PI * t) * std::cos(M_PI * x));

    return 0.0;
}

int main()
{
    unsigned int N = 100,
        T = 1000,
        period = 1000;
    double sigma = 1.0,
        h = 0.01,
        tau = 0.1,
        max = 0,
        gr;

    double *pred = (double*)malloc((N + 1) * sizeof(double));
    double *cur = (double*)malloc((N + 1) * sizeof(double));
    if (pred == nullptr || cur == nullptr){
        free(pred);
        free(cur);
        return -1;
    }

    for (int i = 0; i <= N; i++){
        pred[i] = fi(h * i);
    }
    cur[0] = mu_1(tau);
    cur[N] = mu_2(tau);
    for (unsigned int i = 1; i < N; i++){
        cur[i] = tau * (psi(h * i) + 0.5 * tau / (h * h) * (pred[i - 1] - 2 * pred[i] + pred[i + 1]) + f(h * i, 0.0)) + pred[i];
    }

    if (koleb(N, pred, cur, f, mu_1, mu_2, sigma, h, tau, T, period)){
        return -1;
    }

    //std::cin.get();

    return 0;
}
