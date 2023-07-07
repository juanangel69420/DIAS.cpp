// Thermalising for 2 non-zero coordinates X1_C and X2_C using only imaginary compononets
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <ctime>
#include "eigen/Eigen/Dense"

int start = std::time(nullptr);
const int N = 2;
const int iterations = 1e7;
const int record = 1000;
const double dt = 1e-4;
const double g = 0.05;

typedef std::complex<double> complex;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

matrix commutator(matrix A, matrix B)
{
    return A*B - B*A;
} 

double H(double g, matrix X1_C, matrix X2_C, matrix V1_C, matrix V2_C)
{
    // Compute kinetic energy T
    complex T = 1/(2*g*g) * (V1_C*V1_C + V2_C*V2_C).trace();

    matrix X[2] = {X1_C,X2_C}; 
    matrix commutator_sum;  

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if(i == j)
                continue;
            commutator_sum.noalias() += commutator(X[i],X[j])*(commutator(X[i],X[j])); //can likely be more efficient by less function calls
        }
    }
    complex U = - (1)/(4*g*g) * commutator_sum.trace();
    return std:: abs(T + U);
}

matrix F1_C(matrix X1_C, matrix X2_C)
{
    return commutator(X2_C, commutator(X1_C,X2_C));
}

matrix F2_C(matrix X1_C, matrix X2_C)
{
    return commutator(X1_C, commutator(X2_C,X1_C));
}

matrix gauss_law(
    matrix X1_C, matrix X2_C, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1_C, matrix V2_C, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i ++)
    {
        result.noalias() = commutator(X1_C,V1_C) + commutator(X2_C,V2_C) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}

void update(
    double dt, 
    matrix* X1_C, matrix* X2_C, matrix* X1_C_n, matrix* X2_C_n, 
    matrix* V1_C, matrix* V2_C, matrix* F1_C_0, matrix* F2_C_0, matrix* F1_C_n, matrix* F2_C_n)
{
    *X1_C_n = *X1_C + (*V1_C)*dt + 0.5*(*F1_C_0)*pow(dt,2);
    *X2_C_n = *X2_C + (*V2_C)*dt + 0.5*(*F2_C_0)*pow(dt,2);

    *F1_C_n = F1_C(*X1_C_n,*X2_C_n);
    *F2_C_n = F2_C(*X1_C_n,*X2_C_n);

    *V1_C = *V1_C + 0.5*(*F1_C_0 + *F1_C_n)*dt;
    *V2_C = *V2_C + 0.5*(*F2_C_0 + *F2_C_n)*dt;

    *X1_C = *X1_C_n, *X2_C = *X2_C_n;
    *F1_C_0 = *F1_C_n, *F2_C_0 = *F2_C_n;
}

matrix X1_C_sols[iterations/record], X2_C_sols[iterations/record];
matrix V1_C_sols[iterations/record],V2_C_sols[iterations/record];

int main()
{
    matrix X1_C, X2_C, V1_C, V2_C;
    matrix X1_R, X2_R, V1_R, V2_R;

    std:: mt19937 rng(std::time(nullptr));
    std:: normal_distribution<double> gauss_dist(0,1);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            double z1 = gauss_dist(rng);
            double z2 = gauss_dist(rng);
            X1_C(i,j) = complex(0,z1), X1_C(j,i) = std:: conj(X1_C(i,j));
            X2_C(i,j) = complex(0,z2), X2_C(j,i) = std:: conj(X2_C(i,j));

            X1_R(i,j) = complex(z1,0), X1_R(j,i) = std:: conj(X1_R(i,j));
            X2_R(i,j) = complex(z2,0), X2_R(j,i) = std:: conj(X2_R(i,j));

            if (i == j && i != N-1)
            {
                X1_C(i,j) = complex(z1,0);
                X2_C(i,j) = complex(z2,0);

            }
            if (i == j && i == N-1)
            {
                complex diag_sum_1 = complex(0,0);
                complex diag_sum_2 = complex(0,0);
                complex diag_sum_3 = complex(0,0);
                complex diag_sum_4 = complex(0,0);
                for (int k = 0; k < N-1; k++)
                {
                    diag_sum_1 += X1_C(k,k);
                    diag_sum_2 += X2_C(k,k);
                    diag_sum_3 += X1_R(k,k);
                    diag_sum_4 += X2_R(k,k);
                }
                X1_C(i,j) = - diag_sum_1;
                X2_C(i,j) = - diag_sum_2;
                X1_R(i,j) = - diag_sum_3;
                X2_R(i,j) = - diag_sum_4;
            }
        }
    }

    // Initializing F function at t = 0 for use in the update function
    matrix F1_C_0 = F1_C(X1_C,X2_C);
    matrix F2_C_0 = F2_C(X1_C,X2_C);
    matrix X1_C_n, X2_C_n, F1_C_n, F2_C_n;

    matrix F1_R_0 = F1_C(X1_R,X2_R);
    matrix F2_R_0 = F2_C(X1_R,X2_R);
    matrix X1_R_n, X2_R_n, F1_R_n, F2_R_n;

    for (int i = 0; i < 100000; i++)
    {
        if(i%10000 == 0){
            std:: cout << "X1_C: " << std:: endl << X1_C << std::endl << "X1_R: " << std:: endl << X1_R << std:: endl;
            std:: cout << "X2_C: " << std:: endl << X2_C << std::endl << "X2_R: " << std:: endl << X2_R << std:: endl; 
            std:: cout << "Complex energy: " << H(g,X1_C,X2_C,V1_C,V2_C) << std:: endl;
            std:: cout << "Real Energy: " << H(g,X1_R,X2_R,V1_R,V2_R) << std:: endl;
        }
        update(dt, &X1_C, &X2_C, &X1_C_n, &X2_C_n, &V1_C, &V2_C, &F1_C_0, &F2_C_0, &F1_C_n, &F2_C_n);
        update(dt, &X1_R, &X2_R, &X1_R_n, &X2_R_n, &V1_R, &V2_R, &F1_R_0, &F2_R_0, &F1_R_n, &F2_R_n);

        if (std:: isnan(real(X1_C(0,0)))){
            break;
        }
        /*
        if (i%100000 == 0)
        {
            std:: cout << i << std:: endl;
            std:: cout << H(g,X1_C,X2_C,V1_C,V2_C) << std:: endl;
        }
        if(i%record == 0)
        {
            X1_C_sols[i/record] = X1_C, X2_C_sols[i/record] = X2_C, V1_C_sols[i/record] = V1_C, V2_C_sols[i/record] = V2_C;
        }
        */
    }

    std:: cout << X1_C << std:: endl << X1_C.adjoint();
    std:: fstream X_file("C:/Users/cilli/DIAS.cpp/testing/X_file.txt",std:: ios:: out);
    std:: fstream V_file("C:/Users/cilli/DIAS.cpp/testing/V_file.txt",std:: ios:: out);
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X1_C_sols[k](i,j);
                    V_file << V1_C_sols[k](i,j);
                }
            }
            X_file << ":";
            V_file << ":";
        }
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X2_C_sols[k](i,j);
                    V_file << V2_C_sols[k](i,j);
                }
            }
            X_file << ":";
            V_file << ":";
        }
    return 0;
}