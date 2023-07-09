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
    double x1 = 100;
    double x2 = 200;
    std:: cout << x1/x2;
}