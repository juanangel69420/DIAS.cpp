#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include <fstream>
#include "eigen/Eigen/Dense"

int start = std::time(nullptr);
const int N = 2;
const int iterations = 1e6;
const double dt = 1e-4;
const double g = 0.000001;
const int record = 1000;

typedef std::complex<double> complex;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;

matrix commutator(matrix A, matrix B)
{
    return A*B - B*A;
} 

double H(
    double g, 
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    // Compute kinetic energy T
    complex T = 0.5 * (V1*V1 + V2*V2 + V3*V3 + V4*V4 + V5*V5 + V6*V6 + V7*V7 + V8*V8 + V9*V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum;  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i == j)
                continue;
            commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls
        }
    }
    complex U = - (g*g)/(4) * commutator_sum.trace();
    return std:: abs(T + U);
}

matrix F(int i, matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7,matrix X8, matrix X9)
{
    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9};
    matrix sum;
    for (int k = 0; k < 9; k++)
    {
        if(i-1 == k)
        {
            continue;
        }
        sum += commutator(X[k],commutator(X[i - 1],X[k]));
    }
    return g*g*sum;
}

matrix gauss_law(
    matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7, matrix X8, matrix X9,
    matrix V1, matrix V2, matrix V3, matrix V4, matrix V5, matrix V6, matrix V7, matrix V8, matrix V9)
{
    matrix result;
    for (int i = 0; i < 9; i ++)
    {
        result = commutator(X1,V1) + commutator(X2,V2) + commutator(X3,V3) + commutator(X4,V4) + commutator(X5,V5) +
        commutator(X6,V6) + commutator(X7,V7) + commutator(X8,V8) + commutator(X9,V9);
    }
    return result;
}

void update(
    double dt,
    matrix* X1, matrix* X2, matrix* X3, matrix* X4, matrix* X5, matrix* X6, matrix* X7, matrix* X8, matrix* X9,
    matrix* X1_n, matrix* X2_n, matrix* X3_n, matrix* X4_n, matrix* X5_n, matrix* X6_n, matrix* X7_n, matrix* X8_n, matrix* X9_n,
    matrix* V1, matrix* V2, matrix* V3, matrix* V4, matrix* V5, matrix* V6, matrix* V7, matrix* V8, matrix* V9,
    matrix* F1_0, matrix* F2_0, matrix* F3_0, matrix* F4_0, matrix* F5_0, matrix* F6_0, matrix* F7_0, matrix* F8_0, matrix* F9_0,
    matrix* F1_n, matrix* F2_n, matrix* F3_n, matrix* F4_n, matrix* F5_n, matrix* F6_n, matrix* F7_n, matrix* F8_n, matrix* F9_n)
{
    *X1_n = *X1 + (*V1)*dt + 0.5*(*F1_0)*pow(dt,2);
    *X2_n = *X2 + (*V2)*dt + 0.5*(*F2_0)*pow(dt,2);
    *X3_n = *X3 + (*V3)*dt + 0.5*(*F3_0)*pow(dt,2);
    *X4_n = *X4 + (*V4)*dt + 0.5*(*F4_0)*pow(dt,2);
    *X5_n = *X5 + (*V5)*dt + 0.5*(*F5_0)*pow(dt,2);
    *X6_n = *X6 + (*V6)*dt + 0.5*(*F6_0)*pow(dt,2);
    *X7_n = *X7 + (*V7)*dt + 0.5*(*F7_0)*pow(dt,2);
    *X8_n = *X8 + (*V8)*dt + 0.5*(*F8_0)*pow(dt,2);
    *X9_n = *X9 + (*V9)*dt + 0.5*(*F9_0)*pow(dt,2);

    *F1_n = F(1,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F2_n = F(2,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F3_n = F(3,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F4_n = F(4,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F5_n = F(5,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F6_n = F(6,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F7_n = F(7,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F8_n = F(8,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);
    *F9_n = F(9,*X1_n,*X2_n,*X3_n,*X4_n,*X5_n,*X6_n,*X7_n,*X8_n,*X9_n);

    *V1 = *V1 + 0.5*(*F1_0 + *F1_n)*dt;
    *V2 = *V2 + 0.5*(*F2_0 + *F2_n)*dt;
    *V3 = *V3 + 0.5*(*F3_0 + *F3_n)*dt;
    *V4 = *V4 + 0.5*(*F4_0 + *F4_n)*dt;
    *V5 = *V5 + 0.5*(*F5_0 + *F5_n)*dt;
    *V6 = *V6 + 0.5*(*F6_0 + *F6_n)*dt;
    *V7 = *V7 + 0.5*(*F7_0 + *F7_n)*dt;
    *V8 = *V8 + 0.5*(*F8_0 + *F8_n)*dt;
    *V9 = *V9 + 0.5*(*F9_0 + *F9_n)*dt;

    *X1 = *X1_n, *X2 = *X2_n, *X3 = *X3_n, *X4 = *X4_n, *X5 = *X5_n, *X6 = *X6_n, *X7 = *X7_n, *X8 = *X8_n, *X9 = *X9_n;
    *F1_0 = *F1_n, *F2_0 = *F2_n, *F3_0 = *F3_n, *F4_0 = *F4_n, *F5_0 = *F5_n, *F6_0 = *F6_n, *F7_0 = *F7_n, *F8_0 = *F8_n, *F9_0 = *F9_n;
}

// Declare solution arrays
    matrix X1_sols[record],X2_sols[record],X3_sols[record],X4_sols[record],X5_sols[record],X6_sols[record],X7_sols[record],X8_sols[record],X9_sols[record];
    matrix V1_sols[record],V2_sols[record],V3_sols[record],V4_sols[record],V5_sols[record],V6_sols[record],V7_sols[record],V8_sols[record],V9_sols[record];
    matrix X1_sols_ds[record],X2_sols_ds[record],X3_sols_ds[record],X4_sols_ds[record],X5_sols_ds[record],X6_sols_ds[record],X7_sols_ds[record],X8_sols_ds[record],X9_sols_ds[record];
    matrix V1_sols_ds[record],V2_sols_ds[record],V3_sols_ds[record],V4_sols_ds[record],V5_sols_ds[record],V6_sols_ds[record],V7_sols_ds[record],V8_sols_ds[record],V9_sols_ds[record];

int main()
{
    // Declaring the matrices
    matrix X1, X2, X3, X4, X5, X6, X7, X8, X9, X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds;
    matrix V1, V2, V3, V4, V5, V6, V7, V8, V9, V1_ds, V2_ds, V3_ds, V4_ds, V5_ds, V6_ds, V7_ds, V8_ds, V9_ds;

    // Loading Unperturbed and perturbed coordinates
    matrix X[18];
    matrix X_ds[18];
    std:: fstream unperturbed("Thermalised_branes.txt", std:: ios:: in);
    std:: fstream perturbed("Perturbed_branes.txt", std:: ios:: in);
    for (int i = 0; i < 18; i++)
    {
        for (int k = 0; k < N; k++)
        {
            for (int j = 0; j < N; j++)
            {
                unperturbed >> X[i](k, j);
                perturbed >> X_ds[i](k, j);
            }
        }
    }

    X1 = X[0], X2 = X[1], X3 = X[2], X4 = X[3], X5 = X[4], X6 = X[5], X7 = X[6], X8 = X[7], X9 = X[8];
    V1 = X[9], V2 = X[10], V3 = X[11], V4 = X[12], V5 = X[13], V6 = X[14], V7 = X[15], V8 = X[16], V9 = X[17];
    X1_ds = X_ds[0], X2_ds = X_ds[1], X3_ds = X_ds[2], X4_ds = X_ds[3], X5_ds = X_ds[4], X6 = X_ds[5], X7_ds = X_ds[6], X8_ds = X_ds[7], X9_ds = X_ds[8];
    V1_ds = X_ds[9], V2_ds = X_ds[10], V3_ds = X_ds[11], V4_ds = X_ds[12], V5_ds = X_ds[13], V6_ds = X_ds[14], V7_ds = X_ds[15], V8_ds = X_ds[16], V9_ds = X_ds[17];

    // Initializing F function at t = 0 for use in the update function
    matrix F1_0 = F(1,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F2_0 = F(2,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F3_0 = F(3,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F4_0 = F(4,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F5_0 = F(5,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F6_0 = F(6,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F7_0 = F(7,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F8_0 = F(8,X1,X2,X3,X4,X5,X6,X7,X8,X9);
    matrix F9_0 = F(9,X1,X2,X3,X4,X5,X6,X7,X8,X9);

    matrix F1_0_ds = F(1,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F2_0_ds = F(2,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F3_0_ds = F(3,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F4_0_ds = F(4,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F5_0_ds = F(5,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F6_0_ds = F(6,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F7_0_ds = F(7,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F8_0_ds = F(8,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);
    matrix F9_0_ds = F(9,X1_ds, X2_ds, X3_ds, X4_ds, X5_ds, X6_ds, X7_ds, X8_ds, X9_ds);

    // Initializing X1_n's and F1_n's to pass addresses to the update function
    matrix X1_n, X2_n, X3_n, X4_n, X5_n, X6_n, X7_n, X8_n, X9_n, X1_n_ds, X2_n_ds, X3_n_ds, X4_n_ds, X5_n_ds, X6_n_ds, X7_n_ds, X8_n_ds, X9_n_ds;
    matrix F1_n, F2_n, F3_n, F4_n, F5_n, F6_n, F7_n, F8_n, F9_n, F1_n_ds, F2_n_ds, F3_n_ds, F4_n_ds, F5_n_ds, F6_n_ds, F7_n_ds, F8_n_ds, F9_n_ds;

    // Run update function
    for (int i = 0; i < 1e5; i++)
    {
        update(
            dt,
            &X1,&X2,&X3,&X4,&X5,&X6,&X7,&X8,&X9,
            &X1_n,&X2_n,&X3_n,&X4_n,&X5_n,&X6_n,&X7_n,&X8_n,&X9_n,
            &V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8,&V9,
            &F1_0,&F2_0,&F3_0,&F4_0,&F5_0,&F6_0,&F7_0,&F8_0,&F9_0,
            &F1_n,&F2_n,&F3_n,&F4_n,&F5_n,&F6_n,&F7_n,&F8_n,&F9_n);
    
        update(
            dt,
            &X1_ds,&X2_ds,&X3_ds,&X4_ds,&X5_ds,&X6_ds,&X7_ds,&X8_ds,&X9_ds,
            &X1_n_ds,&X2_n_ds,&X3_n_ds,&X4_n_ds,&X5_n_ds,&X6_n_ds,&X7_n_ds,&X8_n_ds,&X9_n_ds,
            &V1_ds,&V2_ds,&V3_ds,&V4_ds,&V5_ds,&V6_ds,&V7_ds,&V8_ds,&V9_ds,
            &F1_0_ds,&F2_0_ds,&F3_0_ds,&F4_0_ds,&F5_0_ds,&F6_0_ds,&F7_0_ds,&F8_0_ds,&F9_0_ds,
            &F1_n_ds,&F2_n_ds,&F3_n_ds,&F4_n_ds,&F5_n_ds,&F6_n_ds,&F7_n_ds,&F8_n_ds,&F9_n_ds);

        if (i % record == 0)
        {
            X1_sols[i/record] = X1, X2_sols[i/record] = X2, X3_sols[i/record] = X3, X4_sols[i/record] = X4, 
            X5_sols[i/record] = X5, X6_sols[i/record] = X6, X7_sols[i/record] = X7, X8_sols[i/record] = X8, X9_sols[i/record] = X9;
            
            X1_sols_ds[i/record] = X1_ds, X2_sols_ds[i/record] = X2_ds, X3_sols_ds[i/record] = X3_ds, X4_sols_ds[i/record] = X4_ds,
            X5_sols_ds[i/record] = X5_ds, X6_sols_ds[i/record] = X6_ds, X7_sols_ds[i/record] = X7_ds, X8_sols_ds[i/record] = X8_ds,
            X9_sols_ds[i/record] = X9;

            V1_sols[i/record] = V1, V2_sols[i/record] = V2, V3_sols[i/record] = V3, V4_sols[i/record] = V4,
            V5_sols[i/record] = V5, V6_sols[i/record] = V6, V7_sols[i/record] = V7, V8_sols[i/record] = V8,
            V9_sols[i/record] = V9;

            V1_sols_ds[i/record] = V1_ds, V2_sols_ds[i/record] = V2_ds, V3_sols_ds[i/record] = V3_ds, V4_sols_ds[i/record] = V4_ds,
            V5_sols_ds[i/record] = V5_ds, V6_sols_ds[i/record] = V6_ds, V7_sols_ds[i/record] = V7_ds, V8_sols_ds[i/record] = V8_ds,
            V9_sols_ds[i/record] = V9;
        }

        if (i%10000 == 0)
        {
            std:: cout << i << std::endl;
            std:: cout << "time: " << std::time(nullptr) - start << std:: endl;
            std:: cout << H(g,X1,X2,X3,X4,X5,X6,X7,X8,X9,V1,V2,V3,V4,V5,V6,V7,V8,V9) << std:: endl;
            std:: cout << H(g,X1_ds,X2_ds,X3_ds,X4_ds,X5_ds,X6_ds,X7_ds,X8_ds,X9_ds,V1_ds,V2_ds,V3_ds,V4_ds,V5_ds,V6_ds,V7_ds,V8_ds,V9_ds) << std:: endl;
        }
    }
    return 0;
}