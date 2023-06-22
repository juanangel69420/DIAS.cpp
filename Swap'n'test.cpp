#include <iostream>
#include <typeinfo>
#include <cmath>

int main()
{
    double x = sqrt(-2);
    std:: cout << x << std:: endl;
    if(x == nan)
        std:: cout << "x is true";
    return 0;
}