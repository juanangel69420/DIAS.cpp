#include <iostream>
#include <cmath>

int main()
{
    double x = sqrt(-2);
    if(std:: isnan(x))
    {
        std:: cout << "x is nan";
    }
    else
        std:: cout << "x is a number";
    
    return 0;
}