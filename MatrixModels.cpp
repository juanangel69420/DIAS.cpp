#include <iostream>
#include <random>
#include <ctime>
#include <complex>

const int N = 2;

std:: complex<double> (*(createMatrix)())[N]
{
    std:: mt19937 rng(std::time(nullptr));
    std:: normal_distribution<double> gauss_dist(0,1);
    static std::complex<double> arr[N][N];

    for (int i=0; i < N; i++)
    {
        for (int j=0; j < N; j++)
        {
            arr[i][j] = std:: complex<double> (gauss_dist(rng),gauss_dist(rng));
        }
    }
    return arr;
}



int main()
{
    std::complex<double> (*X1)[N], (*X2)[N], (*X3)[N], (*X4)[N], (*X5)[N], (*X6)[N], (*X7)[N], (*X8)[N], (*X9)[N];

    X1 = createMatrix(), X2 = createMatrix(), X3 = createMatrix(), X4 = createMatrix(),
    X5 = createMatrix(), X6 = createMatrix(), X7 = createMatrix(), X8 = createMatrix(), X9 = createMatrix();

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            std:: cout << "Matrix 1: " <<std:: endl << X1[i][j] << std:: endl;
            std:: cout << "Matrix 2: " <<std:: endl << X2[i][j] << std:: endl;
            std:: cout << "Matrix 3: " <<std:: endl << X3[i][j] << std:: endl;
            std:: cout << "Matrix 4: " <<std:: endl << X4[i][j] << std:: endl;
            std:: cout << "Matrix 5: " <<std:: endl << X5[i][j] << std:: endl;
            std:: cout << "Matrix 6: " <<std:: endl << X6[i][j] << std:: endl;
            std:: cout << "Matrix 7: " <<std:: endl << X7[i][j] << std:: endl;
            std:: cout << "Matrix 8: " <<std:: endl << X8[i][j] << std:: endl;
            std:: cout << "Matrix 9: " <<std:: endl << X9[i][j] << std:: endl;
        }
    }
    return 0;
}