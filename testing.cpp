#include <iostream>
#include <string>
#include <cmath>
#include "eigen/Eigen/Dense"

const int N = 2;

typedef std::complex<double> complex;
typedef Eigen:: Matrix<complex, N, N> matrix;
typedef Eigen:: Matrix<complex, 1, N> row;
typedef Eigen:: Matrix<complex, N, 1> col;

int main()
{
    col a = col::Random();
    col b = col::Random();
    std::cout << (a.adjoint()*a).real() << std::endl;
    std::cout << (a.adjoint()*a).trace().real();
}