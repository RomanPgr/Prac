#define _USE_MATH_DEFINES
//#define EIGEN_STACK_ALLOCATION_LIMIT 10000
#include <iostream>
#include <cmath>
#include <cstdio>
#include <string.h>
#include <stdlib.h>
#include "C:\eigen-3.3.8\Eigen\Dense"
#include "C:\eigen-3.3.8\Eigen\SparseCore"

#include "C:\eigen-3.3.8\Eigen\Sparse"
#include "C:\eigen-3.3.8\Eigen\IterativeLinearSolvers"
//#include "C:\eigen-3.3.8\Eigen\SuperLUSupport"
#include "C:\eigen-3.3.8\Eigen\SparseCholesky"




#define M 971
#define M_STAR 30


#define T_CH (1e-8)
#define K_CH (1e-16)
#define N_CH (1e+23)
#define N (M+M_STAR-1)     //Количество уравнений
#define NUMBER_NON_ZERO (N+(4*(M_STAR-2))+3+(3*(M-M_STAR-1))+2+(3*(M_STAR-1)))  //Количество ненулевых элементов в матрице Якоби
#define N1 ((1e+23)/N_CH)  //Начальная концентрация мономера
#define N_G ((1e+25)/N_CH) //Концентрация буферного газа



#define STEP_RUNGE_KUTTA (0.000001)
#define STEP_GIRA (0.01)
#define T_MAX (0.10)
#define CHASTOTA (1)


double k(const unsigned int a, const unsigned int b);
double t_dis(const unsigned int n);
void runge_kutta(Eigen::MatrixXd& result, const Eigen::MatrixXd& begin, const double max = STEP_GIRA);
void f_point(Eigen::MatrixXd& result, const Eigen::MatrixXd& y);
void jacobi(Eigen::MatrixXd& result, const Eigen::MatrixXd& param);
void gir(const Eigen::MatrixXd& begin);

double k(const unsigned int a, const unsigned int b){
	//return something / K_CH
	return 1.0;
}

double t_dis(const unsigned int n){
	double t_dis_2 = 2.1*(1e-14 / T_CH);

	if (n > 2){
		t_dis_2 = std::exp(2.5 * n) * t_dis_2;
	}

	switch (n){
	case 2:
		t_dis_2 = 2.1e-14;
		break;
	case 3:
		t_dis_2 = 1.4e-12;
		break;
	case 4:
		t_dis_2 = 8.9e-11;
		break;
	case 5:
		t_dis_2 = 8.0e-10;
		break;
	case 6:
		t_dis_2 = 7.3e-9;
		break;
	case 7:
		t_dis_2 = 2.2e-8;
		break;
	case 8:
		t_dis_2 = 9.9e-8;
		break;
	case 9:
		t_dis_2 = 3.5e-7;
		break;
	case 10:
		t_dis_2 = 5.9e-7;
		break;
	case 11:
		t_dis_2 = 1.6e-6;
		break;
	case 12:
		t_dis_2 = 2.2e-6;
		break;
	case 13:
		t_dis_2 = 4.7e-6;
		break;
	case 14:
		t_dis_2 = 9.6e-6;
		break;
	case 15:
		t_dis_2 = 1.1e-5;
		break;
	case 16:
		t_dis_2 = 2.0e-5;
		break;
	case 17:
		t_dis_2 = 3.5e-5;
		break;
	case 18:
		t_dis_2 = 3.8e-5;
		break;
	case 19:
		t_dis_2 = 6.0e-5;
		break;
	case 20:
		t_dis_2 = 6.3e-5;
		break;
	case 21:
		t_dis_2 = 9.5e-5;
		break;
	case 22:
		t_dis_2 = 1.4e-4;
		break;
	case 23:
		t_dis_2 = 1.4e-4;
		break;
	case 24:
		t_dis_2 = 2.0e-4;
		break;
	case 25:
		t_dis_2 = 2.7e-4;
		break;
	case 26:
		t_dis_2 = 2.8e-4;
		break;
	case 27:
		t_dis_2 = 3.6e-4;
		break;
	case 28:
		t_dis_2 = 3.5e-4;
		break;
	case 29:
		t_dis_2 = 4.6e-4;
		break;
	case 30:
		t_dis_2 = 5.9e-4;
		break;
	default:
		throw ("n > 30 " + __LINE__);
		break;
	}


	t_dis_2 /= T_CH;
	return t_dis_2;
}

void f_point(Eigen::MatrixXd& result, const Eigen::MatrixXd& y){
	{
		double res = -2 * k(1, 1) * y(0, 0);
		for (unsigned int i = 1; i < M; i++){
			res = res - y(i, 0) * k(i + 1, 1);
		}
		res = res * y(0, 0) + 2 * y(M, 0) / t_dis(2);
		for (unsigned int i = M + 1; i < N; i++){
			res = res + y(i, 0) / t_dis(i - M + 2);
		}
		result(0, 0) = res;
	}

	for (unsigned int i = 1; i < M_STAR - 1; i++){
		result(i, 0) = (k(i + 1, 2) * y(M + i - 1, 0) * N_G) + (y(M + i, 0) / t_dis(i + 2)) - k(i + 1, 1) * y(i, 0) * y(0, 0);
	}
	result(M_STAR - 1, 0) = k(M_STAR, 2) * y(N - 1, 0) * N_G - k(M_STAR, 1) * y(M_STAR - 1, 0) * y(0, 0);

	for (unsigned int i = M_STAR; i < M - 1; i++){
		result(i, 0) = k(i, 1) * y(i - 1, 0) * y(0, 0) - k(i + 1, 1) * y(i, 0) * y(0, 0);
	}
	result(M - 1, 0) = k(M - 1, 1) * y(M - 2, 0) * y(0, 0);

	for (unsigned int i = M; i < N; i++){
		result(i, 0) = k(i - M + 1, 1) * y(i - M, 0) * y(0, 0) - y(i, 0) / t_dis(i - M + 2) - y(i, 0) * k(i - M + 2, 2) * N_G;
	}
}

void jacobi(Eigen::MatrixXd& result, const Eigen::MatrixXd& param){

	{ // для y[0]
		result(0, 0) = -4 * k(1, 1) * param(0, 0);
		for (unsigned int i = 1; i < M; i++){
			result(0, 0) -= k(i + 1, 1) * param(i, 0);
		}

		for (unsigned int i = 1; i < M; i++){
			result(0, i) = -k(i + 1, 1) * param(0, 0);
		}

		result(0, M) = 2 / t_dis(2);
		for (unsigned int i = M + 1; i < N; i++){
			result(0, i) = 1 / t_dis(i - M + 2);
		}
	}

	//для 2 <= n < M_STAR
	for (unsigned int i = 1; i < M_STAR - 1; i++){
		result(i, 0) = -k(i + 1, 1) * param(i, 0);
		for (unsigned int j = 1; j < i; j++){
			result(i, j) = 0.0;
		}
		result(i, i) = -k(i + 1, 1) * param(0, 0);

		for (unsigned int j = i + 1; j < M + i - 1; j++){
			result(i, j) = 0.0;
		}
		result(i, M + i - 1) = k(i + 1, 2) * N_G;
		result(i, M + i - 0) = 1 / t_dis(i + 2);
		for (unsigned int j = M + i + 1; j < N; j++){
			result(i, j) = 0.0;
		}
	}

		{
			//Для n = M_STAR
			for (unsigned int j = 0; j < N; j++){
				result(M_STAR - 1, +j) = 0.0;
			}
			result(M_STAR - 1, 0) = -k(M_STAR, 1) * param(M_STAR - 1, 0);
			result(M_STAR - 1, M_STAR - 1) = -k(M_STAR, 1) * param(0, 0);
			result(M_STAR - 1, N - 1) = k(M_STAR, 2) * N_G;
		}

		//для M_STAR < n < M
		for (unsigned int i = M_STAR; i < M - 1; i++){
			result(i, 0) = k(i, 1) * param(i - 1, 0) - k(i + 1, 1) * param(i, 0);
			for (unsigned int j = 1; j < i - 1; j++){
				result(i, j) = 0.0;
			}
			result(i, i - 1) = k(i, 1) * param(0, 0);
			result(i, i) = -k(i + 1, 1) * param(0, 0);
			for (unsigned int j = i + 1; j < N; j++){
				result(i, j) = 0.0;
			}
		}


		{
			//для n = M
			result(M - 1, 0) = k(M - 1, 1) * param((M - 1) - 1, 0);
			for (unsigned int j = 1; j < M - 1 - 1; j++){
				result(M - 1, j) = 0.0;
			}
			result(M - 1, M - 1 - 1) = k(M - 1, 1) * param(0, 0);
			for (unsigned int j = M - 1; j < N; j++){
				result(M - 1, j) = 0.0;
			}
		}

		{
			//n = 2 для N_star
			result(M, 0) = 2 * k(1, 1) * param(0, 0);
			for (unsigned int j = 1; j < M; j++){
				result(M, j) = 0.0;
			}
			result(M, M) = -1 / t_dis(2) - k(2, 2) * N_G;
			for (unsigned int j = M + 1; j < N; j++){
				result(M, j) = 0.0;
			}
		}

		for (unsigned int i = 3; i <= M_STAR; i++){
			result(i + M - 2, 0) = k(i - 1, 1) * param(i - 2, 0);
			for (unsigned int j = 1; j < i - 2; j++){
				result(i - 2 + M, j) = 0.0;
			}
			result(i - 2 + M, i - 2) = k(i - 1, 1) * param(0, 0);

			for (unsigned int j = i - 1; j < M + i - 2; j++){
				result(i - 2 + M, j) = 0.0;
			}
			result(i - 2 + M, M + i - 2) = -1 / t_dis(i) - k(i, 2) * N_G;
			for (unsigned int j = M + i - 1; j < N; j++){
				result(i - 2 + M, j) = 0.0;
			}
		}
}


void jacobi(Eigen::SparseMatrix<double>& result, const Eigen::MatrixXd& param){

	{ // для y[0]
		result.coeffRef(0, 0) = -4 * k(1, 1) * param(0, 0);
		for (unsigned int i = 1; i < M; i++){
			result.coeffRef(0, 0) -= k(i + 1, 1) * param(i, 0);
		}

		for (unsigned int i = 1; i < M; i++){
			result.coeffRef(0, i) = -k(i + 1, 1) * param(0, 0);
		}

		result.coeffRef(0, M) = 2 / t_dis(2);
		for (unsigned int i = M + 1; i < N; i++){
			result.coeffRef(0, i) = 1 / t_dis(i - M + 2);
		}
	}

	//для 2 <= n < M_STAR
	for (unsigned int i = 1; i < M_STAR - 1; i++){
		result.coeffRef(i, 0) = -k(i + 1, 1) * param(i, 0);
		for (unsigned int j = 1; j < i; j++){
			result.coeffRef(i, j) = 0.0;
		}
		result.coeffRef(i, i) = -k(i + 1, 1) * param(0, 0);

		for (unsigned int j = i + 1; j < M + i - 1; j++){
			result.coeffRef(i, j) = 0.0;
		}
		result.coeffRef(i, M + i - 1) = k(i + 1, 2) * N_G;
		result.coeffRef(i, M + i - 0) = 1 / t_dis(i + 2);
		for (unsigned int j = M + i + 1; j < N; j++){
			result.coeffRef(i, j) = 0.0;
		}
	}

		{
			//Для n = M_STAR
			for (unsigned int j = 0; j < N; j++){
				result.coeffRef(M_STAR - 1, +j) = 0.0;
			}
			result.coeffRef(M_STAR - 1, 0) = -k(M_STAR, 1) * param(M_STAR - 1, 0);
			result.coeffRef(M_STAR - 1, M_STAR - 1) = -k(M_STAR, 1) * param(0, 0);
			result.coeffRef(M_STAR - 1, N - 1) = k(M_STAR, 2) * N_G;
		}

		//для M_STAR < n < M
		for (unsigned int i = M_STAR; i < M - 1; i++){
			result.coeffRef(i, 0) = k(i, 1) * param(i - 1, 0) - k(i + 1, 1) * param(i, 0);
			for (unsigned int j = 1; j < i - 1; j++){
				result.coeffRef(i, j) = 0.0;
			}
			result.coeffRef(i, i - 1) = k(i, 1) * param(0, 0);
			result.coeffRef(i, i) = -k(i + 1, 1) * param(0, 0);
			for (unsigned int j = i + 1; j < N; j++){
				result.coeffRef(i, j) = 0.0;
			}
		}


		{
			//для n = M
			result.coeffRef(M - 1, 0) = k(M - 1, 1) * param((M - 1) - 1, 0);
			for (unsigned int j = 1; j < M - 1 - 1; j++){
				result.coeffRef(M - 1, j) = 0.0;
			}
			result.coeffRef(M - 1, M - 1 - 1) = k(M - 1, 1) * param(0, 0);
			for (unsigned int j = M - 1; j < N; j++){
				result.coeffRef(M - 1, j) = 0.0;
			}
		}

		{
			//n = 2 для N_star
			result.coeffRef(M, 0) = 2 * k(1, 1) * param(0, 0);
			for (unsigned int j = 1; j < M; j++){
				result.coeffRef(M, j) = 0.0;
			}
			result.coeffRef(M, M) = -1 / t_dis(2) - k(2, 2) * N_G;
			for (unsigned int j = M + 1; j < N; j++){
				result.coeffRef(M, j) = 0.0;
			}
		}

		for (unsigned int i = 3; i <= M_STAR; i++){
			result.coeffRef(i + M - 2, 0) = k(i - 1, 1) * param(i - 2, 0);
			for (unsigned int j = 1; j < i - 2; j++){
				result.coeffRef(i - 2 + M, j) = 0.0;
			}
			result.coeffRef(i - 2 + M, i - 2) = k(i - 1, 1) * param(0, 0);

			for (unsigned int j = i - 1; j < M + i - 2; j++){
				result.coeffRef(i - 2 + M, j) = 0.0;
			}
			result.coeffRef(i - 2 + M, M + i - 2) = -1 / t_dis(i) - k(i, 2) * N_G;
			for (unsigned int j = M + i - 1; j < N; j++){
				result.coeffRef(i - 2 + M, j) = 0.0;
			}
		}
}

void runge_kutta(Eigen::MatrixXd& result, const Eigen::MatrixXd& begin, const double max){

	const unsigned int ITER = (unsigned int)(max / STEP_RUNGE_KUTTA);

	Eigen::MatrixXd k1(N, 1);
	Eigen::MatrixXd k2(N, 1);
	Eigen::MatrixXd k3(N, 1);
	Eigen::MatrixXd k4(N, 1);
	Eigen::MatrixXd current(begin);

	for (unsigned int i = 1; i <= ITER; i++){
		f_point(k1, current);
		f_point(k2, current + k1 * (STEP_RUNGE_KUTTA / 2));
		f_point(k3, current + k2 * (STEP_RUNGE_KUTTA / 2));
		f_point(k4, current + k3 * STEP_RUNGE_KUTTA);

		current = current + (STEP_RUNGE_KUTTA / 6) * (k1 + 2 * k2 + 2 * k3 + k4);

		//if (CHASTOTA != 0 && i % CHASTOTA == 0){
		//	for (unsigned int m = 0; m < N; m++){
		//		std::cout << current(m, 0) << ' ';
		//	}
		//	std::cout << std::endl;
		//}
	}

	for (unsigned int i = 0; i < N; i++){
		result(i, 0) = current(i, 0);
	}
}



void gir(const Eigen::MatrixXd& begin){

	const unsigned int ITER = (unsigned int)(T_MAX / STEP_GIRA);

	Eigen::MatrixXd pred(begin);
	Eigen::MatrixXd current(N, 1);
	Eigen::MatrixXd next(N, 1);
	Eigen::MatrixXd f_point1(N, 1);
	Eigen::MatrixXd temp(N, 1);
	Eigen::MatrixXd matrix_jacobi(N, N);
	Eigen::MatrixXd E(N, N);
	for (unsigned int i = 0; i < N; i++){
		for (unsigned int j = 0; j < N; j++){
			E(i, j) = ((i != j) ? 0.0 : 1.0);
		}
	}

	runge_kutta(current, pred);

	if (CHASTOTA == 1){
		for (unsigned int m = 0; m < N; m++){
			std::cout << current(m, 0) << ' ';
		}
		std::cout << std::endl;
	}


	for (unsigned int i = 2; i < ITER; i++){

		f_point(f_point1, current);
		jacobi(matrix_jacobi, current);
		pred = (4 / 3) * current - (1 / 3) * pred + ((STEP_GIRA * 2 / 3) * ((E - STEP_GIRA * matrix_jacobi).inverse() * f_point1));

		//Следующие 3 строки строки ужасны. Нужно просто поменять указатели, но как?
		//memcpy(temp, current.data(), N * sizeof(double));
		//memcpy(current.data(), pred.data(), N * sizeof(double));
		//memcpy(pred.data(), temp, N * sizeof(double));

		temp = current;
		current = pred;
		pred = temp;

		if (CHASTOTA != 0 && i % CHASTOTA == 0){
			for (unsigned int m = 0; m < N; m++){
				std::cout << current(m, 0) << ' ';
			}
			std::cout << std::endl;// '\n';
		}
	}

	std::cout.flush();
}



//void gir(const Eigen::MatrixXd& begin){
//
//	const unsigned int ITER = (unsigned int)(T_MAX / STEP_GIRA);
//
//	Eigen::MatrixXd pred(begin);
//	Eigen::MatrixXd current(N, 1);
//	Eigen::MatrixXd next(N, 1);
//	Eigen::MatrixXd f_point1(N, 1);
//	Eigen::MatrixXd temp(N, 1);
//	Eigen::SparseMatrix<double> matrix_jacobi(N, N);
//	matrix_jacobi.reserve(NUMBER_NON_ZERO);
//	Eigen::SparseMatrix<double> E(N, N);
//	E.setIdentity();
//	//for (unsigned int i = 0; i < N; i++){
//	//	for (unsigned int j = 0; j < N; j++){
//	//		E(i, j) = ((i != j) ? 0.0 : 1.0);
//	//	}
//	//}
//
//	runge_kutta(current, pred);
//
//	if (CHASTOTA == 1){
//		for (unsigned int m = 0; m < N; m++){
//			std::cout << current(m, 0) << ' ';
//		}
//		std::cout << std::endl;
//	}
//
//	//Eigen::VectorXd solve(N);
//
//
//	for (unsigned int i = 2; i < ITER; i++){
//
//		f_point(f_point1, current);
//		jacobi(matrix_jacobi, current);
//		/*matrix_jacobi *= (-STEP_GIRA);
//		for (unsigned int j = 0; j < N; j++){
//			matrix_jacobi.coeffRef(j, j) += 1.0;
//		}*/
//		matrix_jacobi = E - STEP_GIRA * matrix_jacobi;
//		Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
//		solver.compute(matrix_jacobi);
//		f_point1 = solver.solve(f_point1);
//		pred = (4 / 3) * current - (1 / 3) * pred + ((STEP_GIRA * 2 / 3) * f_point1);
//		//pred = (4 / 3) * current - (1 / 3) * pred + ((STEP_GIRA * 2 / 3) * ((/*E - STEP_GIRA * */matrix_jacobi).inverse() * f_point1));
//
//		//Следующие 3 строки строки ужасны. Нужно просто поменять указатели, но как?
//		//memcpy(temp, current.data(), N * sizeof(double));
//		//memcpy(current.data(), pred.data(), N * sizeof(double));
//		//memcpy(pred.data(), temp, N * sizeof(double));
//
//		current.swap(pred);
//		//temp = current;
//		//current = pred;
//		//pred = temp;
//
//		if (CHASTOTA != 0 && i % CHASTOTA == 0){
//			for (unsigned int m = 0; m < N; m++){
//				std::cout << current(m, 0) << ' ';
//			}
//			std::cout << std::endl;// '\n';
//		}
//	}
//
//	std::cout.flush();
//}






int main(){
	Eigen::MatrixXd next(N, 1);
	Eigen::MatrixXd begin(N, 1);
	begin(0, 0) = N1; std::cout << begin(0, 0) << ' ';
	for (unsigned int i = 1; i < N; i++){
		begin(i, 0) = 0.0;
		std::cout << begin(i, 0) << ' ';
	}
	std::cout << '\n';

	//runge_kutta(next, begin, T_MAX);

	gir(begin);

	std::cout.flush();
	//std::cin.get();
	return 0;
}