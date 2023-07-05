#include <iostream>
#include <random>
#include <string>
#include <ctime>
#include <cmath>
#include <fstream>
#include "eigen/Eigen/Dense"

int start = std::time(nullptr);
const int N = 2;
const int iterations = 1e7;
const double dt = 1e-4;
const double g = 0.1;

typedef Eigen:: Matrix<double, N, N> matrix;

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
    double T = 1/(2*pow(g,2)) * (V1*V1 + V2*V2 + V3*V3 + V4*V4 + V5*V5 + V6*V6 + V7*V7 + V8*V8 + V9*V9).trace();

    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9}; 

    matrix commutator_sum = matrix::Zero();  
    for (int i = 0; i < 9; i++)
    {
        for (int j = 0; j < 9; j++)
        {
            if(i == j)
                continue;
            commutator_sum += commutator(X[i],X[j])*commutator(X[i],X[j]); //can likely be more efficient by less function calls
        }
    }
    
    double U = - 1/(4*pow(g,2)) * commutator_sum.trace();
    return T + U;
}

matrix F(int i, matrix X1, matrix X2, matrix X3, matrix X4, matrix X5, matrix X6, matrix X7,matrix X8, matrix X9)
{
    matrix X[9] = {X1,X2,X3,X4,X5,X6,X7,X8,X9};
    matrix sum = matrix::Zero();
    for (int k = 0; k < 9; k++)
    {
        if(i-1 == k)
        {
            continue;
        }
        sum += commutator(X[k],commutator(X[i - 1],X[k]));
    }
    return sum;
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

int main()
{
    // Declaring the matrices
    matrix X1, X2, X3, X4, X5, X6, X7, X8, X9, V1, V2, V3, V4, V5, V6, V7, V8, V9;

    std:: fstream thermalised_coordinates("Thermalised_branes.txt",std:: ios::in);
    std:: string line;
    getline(thermalised_coordinates, line);
    line >> X1;
    

    /*
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

    // Initializing X1_n's and F1_n's to pass addresses to the update function
    matrix X1_n, X2_n, X3_n, X4_n, X5_n, X6_n, X7_n, X8_n, X9_n;
    matrix F1_n, F2_n, F3_n, F4_n, F5_n, F6_n, F7_n, F8_n, F9_n;

    // Run update function
    for (int i = 0; i < iterations; i++)
    {
        update(
            dt,
            &X1,&X2,&X3,&X4,&X5,&X6,&X7,&X8,&X9,
            &X1_n,&X2_n,&X3_n,&X4_n,&X5_n,&X6_n,&X7_n,&X8_n,&X9_n,
            &V1,&V2,&V3,&V4,&V5,&V6,&V7,&V8,&V9,
            &F1_0,&F2_0,&F3_0,&F4_0,&F5_0,&F6_0,&F7_0,&F8_0,&F9_0,
            &F1_n,&F2_n,&F3_n,&F4_n,&F5_n,&F6_n,&F7_n,&F8_n,&F9_n);
        
        if (i%50000 == 0)
        {
            std:: cout << i << std::endl;
            std:: cout << "time: " << std::time(nullptr) - start << std:: endl;
            std:: cout << H(g,X1,X2,X3,X4,X5,X6,X7,X8,X9,V1,V2,V3,V4,V5,V6,V7,V8,V9) << std:: endl;
            std:: cout << X1 << "\n" << V1 << "\n";
        }
    }
    */
    return 0;
}