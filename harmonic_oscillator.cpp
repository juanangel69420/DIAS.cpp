#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "eigen/Eigen/Dense"

typedef std::complex<double> complex;
double dt = 1e-4;
const int iterations = 1e6;

void update(double dt, complex* z, complex* z_n, complex* v)
{
    *z_n = *z + (*v)*dt - 0.5*(*z)*pow(dt,2);
    *v = *v + 0.5*(- *z_n - *z)*dt;
    *z = *z_n;
}

double H(complex z,complex v)
{
    double E = 0.5 * (pow(v.real(),2) + pow(v.imag(),2))  + 0.5 * (pow(z.real(),2) + pow(z.imag(),2));
    return E;
}

complex sols[iterations];

int main()
{

    complex v = complex(0,0);
    complex z = complex(10,15);
    complex z_n = complex(0,0);

    for (int i = 0; i < iterations; i++)
    {
        sols[i] = z;
        update(dt, &z, &z_n, &v);
        if(i%10000 == 0)
        {
            std::cout << z << std::endl;
        }
    }
    std::fstream real_comps("C:/Users/cilli/DIAS.cpp/harmonic_oscillator_real.txt",std::ios::out);
    std::fstream imag_comps("C:/Users/cilli/DIAS.cpp/harmonic_oscillator_imag.txt",std::ios::out);
    for(int i = 0; i < iterations; i++)
    {
        if(real_comps.is_open())
        {
            real_comps << sols[i].real() << ",";
            imag_comps << sols[i].imag() << ",";
        }
    }
    real_comps << 10;
    imag_comps << 10;
    real_comps.close();
    imag_comps.close();
    return 0;
}