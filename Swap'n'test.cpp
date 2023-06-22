#include <iostream>
#include <cmath>

int main()
{
    double a[2][2] = {1,2,3,4};
    for (int i = 0; i < 2;i++)
    {
        for (int j = 0; j < 2; j++)
        {
            std:: cout << a[i][j] << std:: endl;
        }
    }
    return 0;
}