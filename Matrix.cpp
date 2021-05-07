
#include <iostream>
#include "Matrix.h"
#include <string>


Matrix::Matrix() {
    this->data = new double* [1]{ new double {1.0} };
    this->colsCount = 1;
    this->rowsCount = 1;
}

Matrix::Matrix(int Rows, int Cols) {
    this->data = new double* [Rows];
    for (int i = 0; i < Rows; i++) {
        this->data[i] = new double[Cols];
        for (int j = 0; j < Cols; j++) {
            if (i == j) {
                this->data[i][j] = 1.0;
            }
            else {
                this->data[i][j] = 0.0;
            }
        }        
    }
    this->rowsCount = Rows;
    this->colsCount = Cols;
}

Matrix::Matrix(Matrix&& other) {
    

    std::swap(this->data, other.data);
    std::swap(this->colsCount, other.colsCount);
    std::swap(this->rowsCount, other.rowsCount);
}

Matrix::Matrix(Matrix& other) {
    
    
    this->data = new double* [other.rowsCount];
    for (int row = 0; row < other.rowsCount; row++) {
        this->data[row] = new double[other.colsCount];
        for (int col = 0; col < other.colsCount; col++)
            this->data[row][col] = other.data[row][col];
    }
    this->colsCount = other.colsCount;
    this->rowsCount = other.rowsCount;
}


int Matrix::getRowsCount() {
    return this->rowsCount;
}
int Matrix::getColsCount() {
    return this->colsCount;
}

Matrix Matrix::operator+(Matrix& right) {
    if (!(right.getColsCount() == this->colsCount && right.getRowsCount() == this->rowsCount)) {
        throw std::string("operator+ : sizes of two matrices don't match up.");
    }
    else {         
        Matrix result(this->rowsCount, this->colsCount);
        for(int row = 0; row < this->rowsCount; row++) {            
            for (int col = 0; col < this->colsCount; col++) {
                result[row][col] = (*this)[row][col] + right[row][col];
            }           
        }
        return Matrix(result);
    }
}
   
Matrix Matrix::operator*(Matrix& right) {
    if (!(right.getRowsCount() == this->colsCount)) {
        throw std::string("operator*: Cols != Rows");
    }
    else {
        Matrix result = filledMatrix(this->rowsCount, right.getColsCount(),0.0);;
        for (int row = 0; row < this->rowsCount; row++) {            
            for (int col = 0; col < right.getColsCount(); col++) {                
                for (int k = 0; k < right.getRowsCount(); k++)
                    result[row][col] += this->data[row][k] * right[k][col];               
            }
            
        }
        return Matrix(result);
    }
}

Matrix operator* (double left,  Matrix& right) {
    return right * left;
}

Matrix Matrix::operator*(double right) {
    Matrix result(this->rowsCount, this->colsCount);
    for (int row = 0; row < this->rowsCount; row++) {           
        for (int col = 0; col < this->colsCount; col++) {
            result[row][col] *= right;
        }           
    }
    return Matrix(result);
}

Matrix Matrix::operator-(Matrix& right) {
    if (!(right.getColsCount() == this->colsCount && right.getRowsCount() == this->rowsCount)) {
        throw std::string("operator+ : sizes of two matrices don't match up.");
    }
    else {
        Matrix result(this->rowsCount, this->colsCount);
        for (int row = 0; row < this->rowsCount; row++) {
            for (int col = 0; col < this->colsCount; col++) {
                result[row][col] = (*this)[row][col] - right[row][col];
            }
        }
        return Matrix(result);
    }
}

double* Matrix::operator[](int index){
    return this->data[index];
}

Matrix& Matrix::operator=(const Matrix& right) {
    if (this == &right) return *this;

    for (int row = 0; row < this->rowsCount; row++) {
        delete[] this->data[row];
    }
    delete[] this->data;

    this->data = new double* [right.rowsCount];
    for (int row = 0; row < right.rowsCount; row++) {
        this->data[row] = new double [right.colsCount];
        for(int col = 0; col < right.colsCount; col++)
            this->data[row][col] = right.data[row][col];
    }
    this->colsCount = right.colsCount;
    this->rowsCount = right.rowsCount;


    return *this;
}
  
Matrix& Matrix::operator=(Matrix&& right) {
    if (this == &right) return *this;

    std::swap(this->data, right.data);
    std::swap(this->colsCount, right.colsCount);
    std::swap(this->rowsCount, right.rowsCount);

   
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Matrix& M) {
    std::string result = "";
    for (int row = 0; row < M.rowsCount; row++) {
        for (int col = 0; col < M.colsCount; col++) {
            result += "\t" + std::to_string(M.data[row][col]);
        }
        result += "\n";
    }
    os << result;
    return os;
}

double Matrix::norm() {
    double result = 0;
    for (int i = 0; i < this->rowsCount; i++) {            
        for (int j = 0; j < this->colsCount; j++) {
            result += this->data[i][j] * this->data[i][j];
        }
            
    }
    return result;
}

Matrix Matrix::tridiagonalMatrix(double a1, double a2, double a3, int N) {
    Matrix result(N, N);
    for (int i = 0; i < N; i++) {
        std::vector<double> row;
        for (int j = 0; j < N; j++) {
            if (i == j) {
                result[i][j] = a1;
            }
            else if (i == j + 1 || i + 1 == j) {
                result[i][j] = a2;
            }
            else if (i == j + 2 || i + 2 == j) {
                result[i][j] = a3;
            }
            else {
                result[i][j] = 0.0;
            }
        }
        
    }
    return Matrix(result);
}

Matrix Matrix::filledMatrix(int rows, int cols, double value) {
    Matrix result(rows, cols);
    for (int i = 0; i < rows; i++) {        
        for (int j= 0; j < cols; j++) {
            result[i][j] = value;
        }      
    }
    return result;
}

Matrix Matrix::Jacobi_Method(Matrix& A, Matrix& b, int maxIter, double eps, int* iterCount, bool printIterations) {
    if (A.getColsCount() != A.getRowsCount()) {
        throw std::string("Jacobi_Method: matrix 'A' have to be a square matrix");
    }
    Matrix x_prev = Matrix::filledMatrix(b.getRowsCount(), 1, 1.0);
    Matrix x_curr = Matrix::filledMatrix(b.getRowsCount(), 1, 1.0);
    double norm_res = 0;
    bool flag = maxIter > 0;
    int k = 0;
    do{
            
        for (int i = 0; i < A.getRowsCount(); i++) {
            double sum = 0;                
            for (int j = 0; j < A.getColsCount(); j++) {
                if (j != i) {
                    sum += A[i][j] * x_prev[j][0];
                }
            }                
            x_curr[i][0] = (b[i][0] - sum)/A[i][i];

        }
        x_prev = x_curr;
        
        ;
        norm_res =((A * x_curr) - b).norm();
        if (printIterations) {
            std::cout << "Jacobi method: norm_res = " << norm_res << std::endl;
        }
        k++;
    } while (norm_res > eps && ((flag * maxIter > k)||(!flag)));
    *iterCount = k;
    return x_curr;
}

Matrix Matrix::Gauss_Seidel_Method(Matrix& A, Matrix& b, int maxIter, double eps, int* iterCount, bool printIterations) {
    if (A.getColsCount() != A.getRowsCount()) {
        throw std::string("Gauss_Seidel_Method: Matrix 'A' have to be a square matrix");
    }
    Matrix x_prev = Matrix::filledMatrix(b.getRowsCount(), 1, 1.0);
    Matrix x_curr = Matrix::filledMatrix(b.getRowsCount(), 1, 1.0);
    double norm_res = 0;
    bool flag = maxIter > 0;
    int k = 0;
    do {

        for (int i = 0; i < A.getRowsCount(); i++) {
            double sum = 0;
            int j = 0;
            for ( j= 0; j < i; j++) {                   
                sum += A[i][j] * x_curr[j][0];                    
            }
            for ( j = i+1; j < A.getColsCount(); j++) {
                sum += A[i][j] * x_prev[j][0];
            }
            x_curr[i][0] = (b[i][0] - sum) / A[i][i];

        }
        x_prev = x_curr;

        norm_res = (A * x_curr - b).norm();
        if (printIterations) {
            std::cout << "Gauss-Seidel method: norm_res = " << norm_res << std::endl;
        }
        k++;
    } while (norm_res > eps && ((flag * maxIter > k) || (!flag)));
    *iterCount = k;
    return x_curr;
}

Matrix Matrix::LU_factorization(Matrix& A, Matrix& b) {
    if (A.getColsCount() != A.getRowsCount()) {
        throw std::string("LU_factorization: matrix 'A' have to be a square matrix");
    }
    Matrix U = A;
    
    Matrix L = Matrix(A.getRowsCount(), A.getColsCount());
    int m = A.getRowsCount();
    for (int k = 0; k < m - 1;k++) {
        for (int j = k + 1; j < m; j++) {
            L[j][k] = U[j][k] / U[k][k];
                
            for (int s = k; s < m; s++) {
                U[j][s] = U[j][s] - L[j][k] * U[k][s];
                    
            }
        }
    }
        
    Matrix y = Matrix::filledMatrix(m, 1, 0);


    for (int i = 0; i < m; i++) {
           
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j][0];
        }
        y[i][0] = (b[i][0] - sum) / L[i][i];
            
    }
      
   
    Matrix x = Matrix::filledMatrix(m, 1, 0);
    for (int i = m-1; i >= 0; i--) {
           
        double sum = 0.0;
        for (int j = i+1; j < m; j++) {
            sum += U[i][j] * x[j][0];
        }
            
        x[i][0] = (y[i][0] - sum) / U[i][i];
          
    }
    return x;
}

Matrix::~Matrix() {
    if (this->data == nullptr) { return; }

    for (int row = 0; row < this->rowsCount; row++) {
        delete[] this->data[row];
    }
    delete[] this->data;
}