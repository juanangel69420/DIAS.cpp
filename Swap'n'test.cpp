#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "eigen/Eigen/Dense"

const int N = 2;
typedef Eigen:: Matrix<std:: complex<double>, N, N> matrix;
typedef std:: complex<double> complex;

int main()
{
    matrix X1,X2;
    matrix X[2] = {X1,X2};
    std:: fstream test_file("testing.txt", std:: ios :: in);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < N; j ++)
        {
            for (int k =0; k < N; k ++)
            {
                test_file >> X[i](j, k);
            }
        }
    }
    std:: cout << X[1];
    return 0;
}