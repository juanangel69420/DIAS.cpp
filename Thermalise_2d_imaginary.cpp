// Thermalising for 2 non-zero coordinates X1 and X2 using only imaginary compononets
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

double H(double g, matrix X1, matrix X2, matrix V1, matrix V2)
{
    // Compute kinetic energy T
    complex T = 1/(2*g*g) * (V1*V1 + V2*V2).trace();

    matrix X[2] = {X1,X2}; 
    matrix commutator_sum;  

    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            if(i == j)
                continue;
            commutator_sum.noalias() += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls
        }
    }
    complex U = - (1)/(4*g*g) * commutator_sum.trace();
    return std:: abs(T + U);
}

matrix F1(matrix X1, matrix X2)
{
    return commutator(X2, commutator(X1,X2));
}

matrix F2(matrix X1, matrix X2)
{
    return commutator(X1, commutator(X2,X1));
}

matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i ++)
    {
        result.noalias() = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}

void update(
    double dt, 
    matrix* X1, matrix* X2, matrix* X1_n, matrix* X2_n, 
    matrix* V1, matrix* V2, matrix* F1_0, matrix* F2_0, matrix* F1_n, matrix* F2_n)
{
    *X1_n = *X1 + (*V1)*dt + 0.5*(*F1_0)*pow(dt,2);
    *X2_n = *X2 + (*V2)*dt + 0.5*(*F2_0)*pow(dt,2);

    *F1_n = F1(*X1_n,*X2_n);
    *F2_n = F2(*X1_n,*X2_n);

    *V1 = *V1 + 0.5*(*F1_0 + *F1_n)*dt;
    *V2 = *V2 + 0.5*(*F2_0 + *F2_n)*dt;

    *X1 = *X1_n, *X2 = *X2_n;
    *F1_0 = *F1_n, *F2_0 = *F2_n;
}

matrix X1_sols[iterations/record], X2_sols[iterations/record];
matrix V1_sols[iterations/record],V2_sols[iterations/record];

int main()
{
    matrix X1, X2, V1, V2;
    std:: mt19937 rng(std::time(nullptr));
    std:: normal_distribution<double> gauss_dist(0,0.1);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            X1(i,j) = complex(0,gauss_dist(rng)), X1(j,i) = std:: conj(X1(i,j));
            X2(i,j) = complex(0,gauss_dist(rng)), X2(j,i) = std:: conj(X2(i,j));
            if (i == N-1 && i == j)
            {
                complex diag_sum_1 = complex(0,0);
                complex diag_sum_2 = complex(0,0);
                for (int k = 0; k < N-1; k++)
                {
                    diag_sum_1 += X1(k,k);
                    diag_sum_2 += X2(k,k);
                }
                X1(i,j) = - diag_sum_1;
                X2(i,j) = - diag_sum_2;
            }
        }
    }

    // Initializing F function at t = 0 for use in the update function
    matrix F1_0 = F1(X1,X2);
    matrix F2_0 = F2(X1,X2);
    matrix X1_n, X2_n, F1_n, F2_n;

    for (int i = 0; i < iterations; i++)
    {
        update(dt, &X1, &X2, &X1_n, &X2_n, &V1, &V2, &F1_0, &F2_0, &F1_n, &F2_n);
        
        if (std:: isnan(real(X1(0,0)))){
            break;
        }
        if (i%100000 == 0)
        {
            std:: cout << i << std:: endl;
            std:: cout << H(g,X1,X2,V1,V2) << std:: endl;
        }
        if(i%record == 0)
        {
            X1_sols[i/record] = X1, X2_sols[i/record] = X2, V1_sols[i/record] = V1, V2_sols[i/record] = V2;
        }
    }

    std:: fstream X_file("C:/Users/cilli/DIAS.cpp/testing/X_file.txt",std:: ios:: out);
    std:: fstream V_file("C:/Users/cilli/DIAS.cpp/testing/V_file.txt",std:: ios:: out);
    for (int k = 0; k < iterations/record; k ++)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                {
                    X_file << X1_sols[k](i,j);
                    V_file << V1_sols[k](i,j);
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
                    X_file << X2_sols[k](i,j);
                    V_file << V2_sols[k](i,j);
                }
            }
            X_file << ":";
            V_file << ":";
        }
    return 0;
}