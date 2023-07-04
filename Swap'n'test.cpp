#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <complex>
#include <fstream>
#include "eigen/Eigen/Dense"

const int N = 2;
typedef Eigen:: Matrix<double, N, N> matrix;

struct strings{
    std:: string str1;
    std:: string str2;
};

strings split_string(std:: string x) // Adjust number of member variables in split screen according to the number of columns in the matrices
{    
    strings result;
    for (int i = 0; i < x.size(); i++)
    {
        if (isspace(x[i]))
            break;
        else
            result.str1 += x[i];
    }

    for (int i = result.str1.size() + 1; i < x.size(); i ++)
    {
        result.str2 += x[i];
    }
    return result;
}

int main()
{
    double values[18*N*N];

    std:: fstream thermalised_coordinates("Thermalised_branes.txt",std:: ios:: in);
    std:: string line;
    int i = 0;
    while(getline(thermalised_coordinates, line) && i < values.length())
    {
        strings current_strings = split_string(line);
        double col1 = stod(current_strings.str1);
        double col2 = stod(current_strings.str2);

    }

    return 0;
}