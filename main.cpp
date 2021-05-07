
#include <iostream>
#include <vector>
#include <time.h>  
#include "Matrix.h"
#include <chrono>


void example1() {
   
    double a1 = 3;
    double a2, a3;
    a2 = a3 = -1.0;
    int N = 50;
    Matrix A = Matrix::tridiagonalMatrix(a1, a2, a3, N);
    Matrix b = Matrix::filledMatrix(N, 1, 0.0);
    for (int i = 0; i < N; i++) {
        b[i][0] = sin(double((i + 1) * (0 + 1)));
    }
    Matrix res;

    int iter_jacobi = 0;
    Matrix x_jacobi = Matrix::Jacobi_Method(A, b, 30, pow(10, -9), &iter_jacobi,true);
    res = (A * x_jacobi) - b;
    std::cout << "Jacobi method, norm(A*x-b)="<< res.norm() << " iter:"<<iter_jacobi << std::endl << std::endl;

    int iter_gauss = 0;
    Matrix x_gauss = Matrix::Gauss_Seidel_Method(A, b, 30, pow(10, -9), &iter_gauss,true);
    res = (A * x_gauss) - b;
    std::cout <<"Gauss-Seidel method, norm(A*x-b)="<<res.norm() <<" iter:" <<iter_gauss << std::endl << std::endl;

    Matrix x_LU = Matrix::LU_factorization(A, b);
    res = (A * x_LU) - b;
    std::cout << "LU factorization, norm(A*x-b)=" << res.norm() << std::endl << std::endl;
}


void example2() {
    
    double a1, a2, a3;
    a1 = 5.0;
    a2 = a3 = -1;
    std::vector<int> Ntab { 100,200,300,400,500 };  
    for(int N : Ntab){
        std::cout << "----------------------------------------------------------------" << std::endl;
        std::cout << "N = " << N << std::endl;
       
        Matrix A = Matrix::tridiagonalMatrix(a1, a2, a3, N);
        
        Matrix b = Matrix::filledMatrix(N, 1, 0.0);
        for (int j = 0; j < N; j++) {
            b[j][0] = sin(double((j + 1) * (0 + 1)));
        }
        
        auto start = std::chrono::high_resolution_clock::now();
        int iter_jacobi = 0;

        Matrix x_jacobi = Matrix::Jacobi_Method(A, b, -1, pow(10, -9), &iter_jacobi, false);

        auto stop = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "Jacobi Method " <<
            "time="<< elapsed.count() / 1000.0 <<
            "[s] "<<"iterations="<<iter_jacobi<<
            " norm(A*x-b)="<<(A* x_jacobi - b).norm() << std::endl;


        start = std::chrono::high_resolution_clock::now();
        int iter_gauss = 0;

        Matrix x_gauss = Matrix::Gauss_Seidel_Method(A, b,-1, pow(10, -9), &iter_gauss,false);

        stop = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"Gauss-Seidel Method " <<
            "time="<< elapsed.count()/1000.0 <<
            "[s] "<<"iterations=" << iter_gauss <<
            " norm(A*x-b)=" << (A * x_gauss - b).norm() << std::endl;


        start = std::chrono::high_resolution_clock::now();

        Matrix x_LU = Matrix::LU_factorization(A, b);

        stop = std::chrono::high_resolution_clock::now();
        elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "LU factorization " <<
            "time="<< elapsed.count() / 1000.0 << "[s] "<<
            "norm(A*x-b)=" << (A * x_gauss - b).norm() << std::endl;

        std::cout << std::endl;
       
    }
}

int main()
{
    
    try {              
        example1();
        example2();
    }
    catch (std::string caught) {
        std::cout << "Error: " << caught << std::endl;
        return 1;
    }
    return 0;
   
   
}
