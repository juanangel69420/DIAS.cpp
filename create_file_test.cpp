#include <iostream>
#include <fstream>
#include "eigen/Eigen/Dense"

const int N = 2;
typedef Eigen:: Matrix<std::complex<double>, N, N> matrix;
typedef std:: complex<double> complex;

int main()
{
    
    complex z1 = complex(2,4); 
    complex z2 = complex(1,9); 
    complex z3 = complex(2,2); 
    complex z4 = complex(1,1);
    complex z5 = complex(1,7); 
    complex z6 = complex(9,1); 
    complex z7 = complex(7,7); 
    complex z8 = complex(9,14);

    matrix X1,X2;
    X1 << z1,z2,z3,z4;
    X2 << z5,z6,z7,z8;
    
    matrix X[2] = {X1,X2};
    std:: fstream test_file("testing.txt", std:: ios :: out);
    for (int i = 0; i < 2; i++)
    {
        for(int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                test_file << X[i](j,k);
            }
        }
    }
   return 0;
}