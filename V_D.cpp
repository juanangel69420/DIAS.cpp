#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "eigen/Eigen/Dense"

const double g = 1;
const int N = 3;
const double dt = 1e-4;
const int iterations = 1e6;

typedef std::complex<double> complex;
typedef Eigen:: Matrix<complex, N, N> matrix;
typedef Eigen:: Matrix<complex, 1, N> row;
typedef Eigen:: Matrix<complex, N, 1> col;

double colabsqr(col a)
{
    complex b = (a.adjoint()*a);
    return b.real();
}

double rowabsqr(row a)
{
    complex b = (a*a.adjoint());
    return b.real();
}

double modsqr(complex a)
{
    return pow(abs(a),2);
}

std:: mt19937 rng(std::time(nullptr));
std:: normal_distribution<long double> gauss_dist(0, 1);

// Initialise the Z variables only and let U be the \dot{Z} coordinates
complex Z12 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z13 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z21 = complex(gauss_dist(rng),gauss_dist(rng));
complex Z23 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z31 = complex(gauss_dist(rng),gauss_dist(rng)); 
complex Z32 = complex(gauss_dist(rng),gauss_dist(rng)); 

col Z41 = col::Random();
col Z42 = col::Random();
col Z43 = col::Random();
row Z14 = row::Random();
row Z24 = row::Random(); 
row Z34 = row::Random();

// Initialise U
complex U12 = complex(0,0);
complex U13 = complex(0,0);
complex U21 = complex(0,0);
complex U23 = complex(0,0);
complex U31 = complex(0,0);
complex U32 = complex(0,0);

col U41 = col::Zero();
col U42 = col::Zero();
col U43 = col::Zero();
row U14 = row::Zero();
row U24 = row::Zero();
row U34 = row::Zero();

// Initialise c's as complex but they may be real ask Swapno if things don't work
complex c1 = complex(gauss_dist(rng),gauss_dist(rng));
complex c2 = complex(gauss_dist(rng),gauss_dist(rng));
complex c3 = complex(gauss_dist(rng),gauss_dist(rng));
complex c4 = complex(gauss_dist(rng),gauss_dist(rng));

// Define the Energy 
double H()
{
    double E;
    return E;
}

// Define the Force functions

complex FZ12(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42)
{
    return Z12*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z12*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g*));
}

complex FZ21(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z24, complex Z32, col Z42)
{
    return - Z21*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z21*(modsqr(Z21) + modsqr(Z23) + rowabsqr(Z24) - modsqr(Z12) - modsqr(Z32) - colabsqr(Z42) - c2/(g*g*));
}

complex FZ13(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43)
{
    return Z13*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    - Z13*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g)); 
}

complex FZ31(
    complex Z12, complex Z13, row Z14, complex Z21, complex Z31, col Z41,
    complex Z23, row Z34, complex Z32, col Z43)
{
    return - Z31*(modsqr(Z12) + modsqr(Z13) + rowabsqr(Z14) - modsqr(Z21) - modsqr(Z31) - colabsqr(Z41) - c1/(g*g))
    + Z31*(modsqr(Z31) + modsqr(Z32) + rowabsqr(Z34) - modsqr(13) - modsqr(Z23) - colabsqr(Z43) - c3/(g*g));
}

complex FZ23(
    complex Z23, complex Z21, row Z24, complex Z12, complex Z32, col Z42,
    complex Z31, complex Z34, complex Z13, col Z43)
{
    Z23*(modsqr(Z23))
}


// Define the Update function
/*
void update(
    double dt, complex* Z12, complex* Z12_n, complex* Z21, complex* Z21_n, complex* Z13, complex* Z13_n, complex* Z31, complex* Z31_n,
    complex* Z23, complex* Z23_n, complex* Z32, complex* Z32_n, col* Z41, col* Z41_n, col* Z42, col* Z42_n, col* Z43,
    col* Z43_n, row* Z14, row* Z14_n, row* Z24, row* Z24_n, row* Z34, row* Z34_n)
{
    *Z12_n = *Z12 + U12*dt + 0.5*(Z12*()-Z12)*pow(dt,2);
}
*/

int main()
{
    return 0;
}