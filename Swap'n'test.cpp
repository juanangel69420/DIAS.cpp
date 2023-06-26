#include <iostream>
#include <random>
#include <ctime>
#include <complex>
typedef std:: complex<double> complex;

const int N = 2;

struct matrix{
    
    int rows;
    int columns;
    complex data[N][N];

    // Default constructor is a random NxN matrix
    matrix() : rows(N), columns(N) 
    {
        // Initialize a random number generator engine and initialize a gaussian distribution with mean 0 and std 1
        std:: mt19937 rng(std::time(nullptr));
        std:: normal_distribution<double> gauss_dist(0,1);

        // Assign random complex numbers to the matrix (Ensure that the matrix is hermitian and traceless)
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                complex z1 = complex(gauss_dist(rng),gauss_dist(rng));
                data[i][j] = z1;
                data[j][i] = std:: conj(z1);
                if (i == (N-1) && j == (N-1))
                {
                    complex current_sum = 0;
                    for (int k = 0; k < N-1; k++)
                    {
                        current_sum += data[k][k];
                    }
                    data[i][j] = -current_sum;
                }
            }
        }
    }
};

int main()
{
    matrix A;
    matrix B;
    for (int i = 0; i < N; i++)
    {
        std:: cout << std::endl;
        for (int j = 0; j < N; j++)
        {
            std:: cout << A.data[i][j];
        }
    }
    std::cout << std:: endl;
    for (int i = 0; i < N; i++)
    {
        std:: cout << std::endl;
        for (int j = 0; j < N; j++)
        {
            std:: cout << B.data[i][j];
        }
    }
}