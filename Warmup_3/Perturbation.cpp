#include <iostream>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cmath>
#include <fstream>

// defining the Hamiltonian
double H_compute(double lam1, double lam2, double x1, double x2, double y1, double y2,double px1, double px2, double py1, double py2)
{
    return 0.5*pow(lam1,2)*pow(x2,6) + 1.5*pow(lam1,2)*pow(x2,4)*pow(y2,2) + 1.5*pow(lam1,2)*pow(x2,2)*pow(y2,4) 
    + 0.5*pow(lam1,2)*pow(y2,6) + lam1*x1*pow(x2,3) - 3*lam1*x1*x2*pow(y2,2) + 3*lam1*pow(x2,2)*y1*y2 - lam1*y1*pow(y2,3)
    +0.5*pow(lam2,2)*pow(x1,6) + 1.5*pow(lam2,2)*pow(x1,4)*pow(y1,2) + 1.5*pow(lam2,2)*pow(x1,2)*pow(y1,4) + 0.5*pow(lam2,2)*pow(y1,6)
    +lam2*pow(x1,3)*x2 + 3*lam2*pow(x1,2)*y1*y2 - 3*lam2*x1*x2*pow(y1,2) -lam2*pow(y1,3)*y2 + 0.5*pow(px1,2) + 0.5*pow(px2,2) 
    +0.5*pow(py1,2) + 0.5*pow(py2,2) + 0.5*pow(x1,2) + 0.5*pow(x2,2) + 0.5*pow(y1,2) + 0.5*pow(y2,2); 
}

double px1_ds(double lam1, double lam2, double H, double x1, double x2, double y1, double y2, double px2, double py1, double py2)
{
    return sqrt(2*(H - 0.5*pow(lam1,2)*pow(x2,6) - 1.5*pow(lam1,2)*pow(x2,4)*pow(y2,2) - 1.5*pow(lam1,2)*pow(x2,2)*pow(y2,4)
    - 0.5*pow(lam1,2)*pow(y2,6) - lam1*x1*pow(x2,3) + 3*lam1*x1*x2*pow(y2,2) -3*lam1*pow(x2,2)*y1*y2 + lam1*y1*pow(y2,3)
    - 0.5*pow(lam2,2)*pow(x1,6) - 1.5*pow(lam2,2)*pow(x1,4)*pow(y1,2) - 1.5*pow(lam2,2)*pow(x1,2)*pow(y1,4) - 0.5*pow(lam2,2)*pow(y1,6)
    - lam2*pow(x1,3)*x2 - 3*lam2*pow(x1,2)*y1*y2 + 3*lam2*x1*x2*pow(y1,2) + lam2*pow(y1,3)*y2 - 0.5*pow(px2,2)
    - 0.5*pow(py1,2) - 0.5*pow(py2,2) - 0.5*pow(x1,2) - 0.5*pow(x2,2) - 0.5*pow(y1,2) - 0.5*pow(y2,2)));
}

int main()
{
    // initialise variables
    double lam1 = 0.0001, lam2 = 0.0002;
    double x1, x2, y1, y2, px1, px2, py1, py2, H;

    // load in current thermalised variables
    std:: fstream initial_conditions("new_coordinates.txt", std:: ios:: in);
    initial_conditions >> x1 >> x2 >> y1 >> y2 >> px1 >> px2 >> py1 >> py2;
    initial_conditions.close();

    // perturb the coordinates
    x1 += 0.1, x2 += 0.1, y1 += 0.1, y2 += 0.1, px2 += 0.1, py1 += 0.1, py2 += 0.1;

    // load in the hamiltonian
    std:: fstream Hamiltonian("CurrentHamiltonian.txt", std:: ios:: in);
    Hamiltonian >> H;
    Hamiltonian.close();
    
    px1 = px1_ds(lam1,lam2,H,x1,x2,y1,y2,px2,py1,py2);

    std:: fstream perturbed("perturbed_coordinates.txt", std:: ios:: out);
    perturbed << x1 << "\n" << x2 << "\n" << y1 << "\n" << y2 << "\n" 
    << px1 << "\n" << px2 << "\n" << py1 << "\n" << py2;
    perturbed.close();

    std:: cout << H_compute(lam1,lam2,x1,x2,y1,y2,px1,px2,py1,py2);
}