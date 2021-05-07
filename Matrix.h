#pragma once
#include <vector>
#include <iostream>
#include <string>



class Matrix {
private:
    double** data = nullptr;
    int rowsCount = 0;
    int colsCount = 0;
public:
    Matrix();
    Matrix(int Rows, int Cols);
    Matrix(Matrix&& other);
    Matrix(Matrix& other);

    int getRowsCount();
    int getColsCount();
    double norm();

    Matrix operator+ (Matrix& right);    
    Matrix operator* (Matrix& right);
    Matrix operator* (double right);
    Matrix operator- (Matrix& right);
    Matrix& operator= (Matrix&& right);
    Matrix& operator= (const Matrix& right);
    double* operator[](int index);

    friend std::ostream& operator<< (std::ostream& os, const Matrix& matrix);
    friend Matrix operator* (double left, const Matrix& right);


    static Matrix filledMatrix(int rows, int cols,double value);    
    static Matrix tridiagonalMatrix(double a1, double a2, double a3, int N);
    static Matrix Jacobi_Method(Matrix& A, Matrix& b, int maxIter, double eps, int * iterCount,bool printIterations);
    static Matrix Gauss_Seidel_Method(Matrix& A, Matrix& b, int maxIter, double eps, int* iterCount, bool printIterations);
    static Matrix LU_factorization(Matrix& A, Matrix& x);

    ~Matrix();

};
