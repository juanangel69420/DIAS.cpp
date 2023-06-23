#include <iostream>
#include <random>
#include <ctime>
#include <complex>

const int N = 2;
std:: complex<double> (*(add)(std::complex<double> A[N][N],std::complex<double> B[N][N]))[N]
{
    static std::complex<double> C[N][N];
    for (int i = 0; i < N; i++)
    {
        
        for (int j = 0; j < N; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }

    }
    return C;
}

std:: complex<double> (*(multiply(std::complex<double> A[N][N],std::complex<double> B[N][N])))[N]
{
    std:: complex<double> C[N][N] = {};
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                C[i][j] += A[i][k] + B[k][i];
            }
        }
    }
    return C;
}

std::complex<double> X1[N][N], X2[N][N], X3[N][N], X4[N][N], X5[N][N], X6[N][N], X7[N][N], X8[N][N], X9[N][N];
int main()
{
    // Create the matrices
    std:: mt19937 rng(std::time(nullptr));
    std:: normal_distribution<double> gauss_dist(0,1); 
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std::complex<double> z1(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z2(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z3(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z4(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z5(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z6(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z7(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z8(gauss_dist(rng),gauss_dist(rng));
            std::complex<double> z9(gauss_dist(rng),gauss_dist(rng));
            X1[i][j] = z1;
            X1[j][i] = std:: conj(z1);
            X2[i][j] = z2;
            X2[j][i] = std:: conj(z2);
            X3[i][j] = z3;
            X3[j][i] = std:: conj(z3);
            X4[i][j] = z4;
            X4[j][i] = std:: conj(z4);
            X5[i][j] = z5;
            X5[j][i] = std:: conj(z5);
            X6[i][j] = z6;
            X6[j][i] = std:: conj(z6);
            X7[i][j] = z7;
            X7[j][i] = std:: conj(z7);
            X8[i][j] = z8;
            X8[j][i] = std:: conj(z8);
            X9[i][j] = z9;
            X9[j][i] = std:: conj(z9);
            // Make sure you're using the && correctly
            if (i == N - 1 && i == j)
            {
                std:: complex<double> current_sum1(0,0);
                std:: complex<double> current_sum2(0,0);
                std:: complex<double> current_sum3(0,0);
                std:: complex<double> current_sum4(0,0);
                std:: complex<double> current_sum5(0,0);
                std:: complex<double> current_sum6(0,0);
                std:: complex<double> current_sum7(0,0);
                std:: complex<double> current_sum8(0,0);
                std:: complex<double> current_sum9(0,0);
                for (int k = 0; k < N - 1; k++)
                {
                    current_sum1 += X1[k][k];
                    current_sum2 += X2[k][k];
                    current_sum3 += X3[k][k];
                    current_sum4 += X4[k][k];
                    current_sum5 += X5[k][k];
                    current_sum6 += X6[k][k];
                    current_sum7 += X7[k][k];
                    current_sum8 += X8[k][k];
                    current_sum9 += X9[k][k];
                }
                X1[i][j] = -current_sum1;
                X2[i][j] = -current_sum2;
                X3[i][j] = -current_sum3;
                X4[i][j] = -current_sum4;
                X5[i][j] = -current_sum5;
                X6[i][j] = -current_sum6;
                X7[i][j] = -current_sum7;
                X8[i][j] = -current_sum8;
                X9[i][j] = -current_sum9;
            }
            
        }
    }

    for (int x = 0; x < N; x++)
    {
        for(int y = 0; y < N; y++)
        {
        std:: cout << X1[x][y] << std:: endl;
        }
    }
    
    std:: cout << "Trace: " << X1[0][0] + X1[1][1];
    return 0;
}