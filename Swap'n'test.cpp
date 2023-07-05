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
    X1 << complex(1,3), complex(3,3), complex(1,9), complex (9,1);
    X2 << complex(5,4), complex(2,1), complex(8,5), complex (5,8);
    matrix X[2] = {X1,X2};
    std:: cout << X[0] <<std:: endl << X[1] << std:: endl;
    X1(1,1) = 1000;
    std::cout <<std:: endl<< "New X1: " << X1 << std:: endl; 
    return 0;
}