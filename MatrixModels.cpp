#include <iostream>
#include <random>
#include <ctime>
#include <complex>


typedef std::complex<double> complex;
const int N = 2;

complex (*(add)(complex A[N][N], complex B[N][N]))[N]
{
    static complex C[N][N];
    for (int i = 0; i < N; i++)
    {
        
        for (int j = 0; j < N; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }

    }
    return C;
}

complex (
    *(add_8)(complex x1[N][N], complex x2[N][N], complex x3[N][N], complex x4[N][N], complex x5[N][N],
        complex x6[N][N], complex x7[N][N], complex x8[N][N]))[N]
{
    static complex C[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i][j] = x1[i][j] + x2[i][j] + x3[i][j] + x4[i][j] + x5[i][j] + x6[i][j] + x7[i][j] + x8[i][j];
        }
    }
    return C;
}

complex (*(subtract)(complex A[N][N], complex B[N][N]))[N]
{
    static complex C[N][N];
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

complex (*(multiply)(complex A[N][N], complex B[N][N]))[N]
{
    static complex C[N][N] = {};
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                C[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
    return C;
}

complex Trace(complex A[N][N])

{
    static complex C = 0;
    for (int i = 0; i < N; i++)
    {
        C += A[i][i];
    }
    
    return C;
}

complex (*(commutator)(complex A[N][N], complex B[N][N]))[N]
{    
    complex (*AB)[N] = multiply(A, B);
    std::cout << "AB:";
    for (int i = 0; i < N; i++)
    {
        std:: cout << std::endl; 
        for (int j = 0; j< N; j++)
        {
            std:: cout << AB[i][j] << " ";
        }
    }
    std::cout<<std::endl;
    std::cout << "B:";
    for (int i = 0; i < N; i++)
    {
        std:: cout << std::endl; 
        for (int j = 0; j< N; j++)
        {
            std:: cout << B[i][j] << " ";
        }
    }
    std:: cout << std:: endl;
    complex (*BA)[N] = multiply(A, B);// This needs to be changed back
    std::cout << "BA:";
    for (int i = 0; i < N; i++)
    {
        std:: cout << std::endl; 
        for (int j = 0; j< N; j++)
        {
            std:: cout << BA[i][j] << " ";
        }
    }
    std:: cout << std:: endl;
    static complex (*C)[N] = subtract(AB,BA);
    std:: cout << "C:" << std:: endl;
    std:: cout << C;
    return C;
}

complex (*(scalar_multiply)(complex A[N][N], complex z))[N]
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N;j++)
        {
            A[i][j] = A[i][j] * z;
        }
    }
    return A;
}

complex (
    *(F)(complex Xi[N][N], complex Xj1[N][N], complex Xj2[N][N], complex Xj3[N][N], complex Xj4[N][N],
        complex Xj5[N][N], complex Xj6[N][N], complex Xj7[N][N], complex Xj8[N][N]))[N]
{
    static complex (*C)[N] = add_8(
        commutator(Xj1,commutator(Xi,Xj1)),commutator(Xj2,commutator(Xi,Xj2)),commutator(Xj3,commutator(Xi,Xj3)),
        commutator(Xj4,commutator(Xi,Xj4)),commutator(Xj5,commutator(Xi,Xj5)),commutator(Xj6,commutator(Xi,Xj6)),
        commutator(Xj7,commutator(Xi,Xj7)),commutator(Xi,Xj8));
    return C;
}

complex X1[N][N], X2[N][N], X3[N][N], X4[N][N], X5[N][N], X6[N][N], X7[N][N], X8[N][N], X9[N][N];

int main()
{
    // Create the matrices
    std:: mt19937 rng(std::time(nullptr));
    std:: normal_distribution<double> gauss_dist(0,1); 
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            complex z1(gauss_dist(rng),gauss_dist(rng));
            complex z2(gauss_dist(rng),gauss_dist(rng));
            complex z3(gauss_dist(rng),gauss_dist(rng));
            complex z4(gauss_dist(rng),gauss_dist(rng));
            complex z5(gauss_dist(rng),gauss_dist(rng));
            complex z6(gauss_dist(rng),gauss_dist(rng));
            complex z7(gauss_dist(rng),gauss_dist(rng));
            complex z8(gauss_dist(rng),gauss_dist(rng));
            complex z9(gauss_dist(rng),gauss_dist(rng));
            X1[i][j] = z1;
            X1[j][i] = std:: conj(z1);
            X2[i][j] = z2;
            X2[j][i] = std:: conj(z2);
            X3[i][j] = z3;
            X3[j][i] = std:: conj(z3);
            X4[i][j] = z4;
            X4[j][i] = std:: conj(z4);
            X5[i][j] = z5;
            X5[j][i] = std:: conj(z5);
            X6[i][j] = z6;
            X6[j][i] = std:: conj(z6);
            X7[i][j] = z7;
            X7[j][i] = std:: conj(z7);
            X8[i][j] = z8;
            X8[j][i] = std:: conj(z8);
            X9[i][j] = z9;
            X9[j][i] = std:: conj(z9);
            // Make sure you're using the && correctly
            if (i == N - 1 && i == j)
            {
                complex current_sum1(0,0);
                complex current_sum2(0,0);
                complex current_sum3(0,0);
                complex current_sum4(0,0);
                complex current_sum5(0,0);
                complex current_sum6(0,0);
                complex current_sum7(0,0);
                complex current_sum8(0,0);
                complex current_sum9(0,0);
                for (int k = 0; k < N - 1; k++)
                {
                    current_sum1 += X1[k][k];
                    current_sum2 += X2[k][k];
                    current_sum3 += X3[k][k];
                    current_sum4 += X4[k][k];
                    current_sum5 += X5[k][k];
                    current_sum6 += X6[k][k];
                    current_sum7 += X7[k][k];
                    current_sum8 += X8[k][k];
                    current_sum9 += X9[k][k];
                }
                X1[i][j] = -current_sum1;
                X2[i][j] = -current_sum2;
                X3[i][j] = -current_sum3;
                X4[i][j] = -current_sum4;
                X5[i][j] = -current_sum5;
                X6[i][j] = -current_sum6;
                X7[i][j] = -current_sum7;
                X8[i][j] = -current_sum8;
                X9[i][j] = -current_sum9;
            }
            
        }
    }

    /*
    X1 = np.array([[1,2],[3,4]])
    X2 = np.array([[2,3],[7,8]])
    X3 = np.array([[4,3],[7,7]])
    X4 = np.array([[8,2],[9,1]])
    X5 = np.array([[1,9],[3,3]])
    X6 = np.array([[1,1],[9,5]])
    X7 = np.array([[5,6],[5,1]])
    X8 = np.array([[3,2],[7,4]])
    X9 = np.array([[2,1],[5,5]])
    */

    /*
    complex X1[N][N] = {{1,2},{3,4}}, X2[N][N] = {{2,3},{7,8}}, X3[N][N] = {{4,3},{7,7}}, X4[N][N] = {{8,2},{9,1}},
    X5[N][N] = {{1,9},{3,3}}, X6[N][N] = {{1,1},{9,5}}, X7[N][N] = {{5,6},{5,1}}, X8[N][N] = {{3,2},{7,4}}, X9[N][N] = {{2,1},{5,5}};

    complex (*C)[N] = F(X1,X2,X3,X4,X5,X6,X7,X8,X9);
    */
    complex A[N][N] = {{1,2},{3,4}};
    complex B[N][N] = {{1,2},{3,5}};
    complex C[N][N] = {{1,1},{1,1}};
    complex (*D)[N] = multiply(A,B);
    complex (*E)[N] = subtract(D,A);
    //complex (*C)[N] = commutator(A,B);
    for (int i = 0; i < N; i++)
    {
        std:: cout << std::endl; 
        for (int j = 0; j< N; j++)
        {
            std:: cout << D[i][j] << " ";
        }
    }
    return 0;
}