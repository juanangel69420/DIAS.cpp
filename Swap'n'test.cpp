#include <iostream>
#include <random>
#include <ctime>
#include <complex>

const int N = 2;

/*
std::complex<double>* add(std::complex<double> A[N][N], std::complex<double> B[N][N])
{
    double* C[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            *C[i][j] = A[i][j] + B[i][j];
        }
    }
    
    return C;
}
*/

double (*(add)(double A[N][N], double B[N][N]))[N]
{
    static double C[N][N];
    for (int i = 0; i < N; i++)
    {
        
        for (int j = 0; j < N; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }

    }
    
    return C;
}

double (*(multiply)(double A[N][N], double B[N][N]))[N]
{
    static double C[N][N] = {};
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                C[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
    return C;
}

int main()
{
    double A[2][2] = {{20,3},{8,9}};
    double B[2][2] = {{1,0},{0,1}};
    double (*C)[N] = multiply(A,B);
    for (int i = 0; i < N; i++)
    {
        std:: cout << std:: endl;
        for (int j = 0; j < N; j++)
        {
            std::cout << C[i][j] << " ";
        }
    }
    return 0;
}