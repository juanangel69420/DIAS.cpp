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
    matrix x1, x2;
    x1 << 1,2,3,4;
    x2 << 1,1,1,1;
    std::cout << x1 + complex(0,1)*x2;
    return 0;
}