#include<iostream>
#include<fstream>

void create_files(int example_number)
{
    for (int i = 1; i <= example_number; i++)
    {
        std:: string path = "C:/Users/cilli/VS programmes/loop_file_testing/run_" + std:: to_string(i) + ".txt";
        std:: string path_ds = "C:/Users/cilli/VS programmes/loop_file_testing/run_ds" + std:: to_string(i) + ".txt"; 
        std:: fstream current_file(path,std::ios::out);
        std:: fstream current_file_ds(path_ds,std::ios::out);
        current_file.close();
        current_file_ds.close();
    }
}

const int iterations = 3e6;
const int record = 1000;
double x1_sol[iterations / record], x2_sol[iterations / record], y1_sol[iterations / record], y2_sol[iterations / record],
px1_sol[iterations / record], px2_sol[iterations / record], py1_sol[iterations / record], py2_sol[iterations / record] = {};

int main()
{    
    // initializing lambdas, time step and iterations
    const double lam1 = 0.0001, lam2 = 0.0002;
    const double dt = 1e-4;

    // initialize initial conditions for coordinates using the thermalised state
    double x1, x2, y1, y2, px1, px2, py1, py2;
    std:: fstream initial_conditions("thermalised_coordinates.txt", std:: ios:: in);
    initial_conditions >> x1 >> x2 >> y1 >> y2 >> px1 >> px2 >> py1 >> py2;
    initial_conditions.close();

    // Create the files needed for storing arrays
    create_files(20);
    const int array_size = 5;
    std:: fstream file_array[array_size];

    for (int i=0; i<array_size/*In future you may need to divide array size by 2 since files only indexed up to 20*/; i++)
    {
        std:: string path = "C:/Users/cilli/VS programmes/Run_files1.0/run_" + std::to_string(i+1) + ".txt";
        file_array[i].open(path, std::ios::out);
        file_array[i] << 1 << "\n" << 2 << "\n" << 3;
        file_array[i].close();
    }

    if (!file_array[0].is_open())
        std:: cout << "File 1 is closed";
    return 0;
}