#include "matrix_class.cpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>


std::vector<double> Matrix::Gauss() {
    if (this->cols != this->rows + 1) {
        throw std::invalid_argument("Wrong dimensions");
    }

    std::vector<std::vector<double>> matrixClone = this->data;
    int n = this->rows;
    int m = this->cols - 1;

    for (int k = 0; k < n; k++) {

        int mainRowIndex = k;
        double mainElement = matrixClone[k][k];

        for (int i = k + 1; i < n; i++) {
            if (std::abs(matrixClone[i][k]) > std::abs(mainElement)) {
                mainElement = matrixClone[i][k];
                mainRowIndex = i;
            }
        }

        if (mainElement == 0.0) {
            throw std::invalid_argument("Division by 0. Unable to solve the system.");
        }

        std::swap(matrixClone[k], matrixClone[mainRowIndex]);

        // Нормалізація поточного рядка
        for (int j = k + 1; j <= m; j++) {
            matrixClone[k][j] /= mainElement;
        }
        // Обнулення нижнього трикутника
        for (int i = k + 1; i < n; i++) {
            double factor = matrixClone[i][k];
            for (int j = k + 1; j <= m; j++) {
                matrixClone[i][j] -= factor * matrixClone[k][j];
            }
        }
    }

    // Зворотній хід (обчислення розв'язку)
    std::vector<double> solution(n, 0.0);

    for (int i = n - 1; i >= 0; i--) {
        solution[i] = matrixClone[i][m];

        for (int j = i + 1; j < n; j++) {
            solution[i] -= matrixClone[i][j] * solution[j];
        }
    }

    return solution;
}


double* Jacobi(Matrix& a, std::vector<double>& y) {
    int n = a.getRows();
    std::vector<double> res(n);

    for (int i = 0; i < n; i++) {
        res[i] = y[i] / a[i][i];
    }

    double eps = 0.0001;
    std::vector<double> Xn(n);

    do {
        for (int i = 0; i < n; i++) {
            Xn[i] = y[i] / a[i][i];
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                } else {
                    Xn[i] -= a[i][j] / a[i][i] * res[j];
                }
            }
        }

        bool flag = true;
        for (int i = 0; i < n; i++) {
            if (std::abs(Xn[i] - res[i]) > eps) {
                flag = false;
                break;
            }
        }

        if (flag) {
            break;
        }

        res = Xn;

    } while (true);

    double* result = new double[n];
    for (int i = 0; i < n; ++i) {
        result[i] = res[i];
    }

    return result;
}

void linearRegression(int n, Matrix& x, Matrix& y) {
    double a, b;

    double xsum = 0, x2sum = 0, ysum = 0, xysum = 0;
    for (int i = 0; i < n; i++) {
        xsum += x[i][0];
        ysum += y[i][0];
        x2sum += x[i][0] * x[i][0];
        xysum += x[i][0] * y[i][0];
    }

    // Обчислення коефіцієнтів a та b
    a = (n * xysum - xsum * ysum) / (n * x2sum - xsum * xsum);
    b = (x2sum * ysum - xsum * xysum) / (x2sum * n - xsum * xsum);

    // Виведення результатів
    std::cout << "S.no" << "\tx" << "\ty(observed)" << "\ty(fitted)" << std::endl;
    std::cout << "---------------------------------" << std::endl;

    for (int i = 0; i < n; i++) {
        double y_fit = a * x[i][0] + b;
        std::cout << i + 1 << ".\t" << x[i][0] << "\t" << y[i][0] << "\t" << y_fit << std::endl;
    }

    std::cout << "\nThe linear fit line is of the form:\n\n" << "y = " << a << "x + " << b << std::endl;
}



int main() {

    Matrix A(3, 4);
    A[0][0] = 2;
    A[0][1] = 1;
    A[0][2] = -1;
    A[1][0] = -3;
    A[1][1] = -1;
    A[1][2] = 2;
    A[2][0] = -2;
    A[2][1] = 1;
    A[2][2] = 2;
    A[0][3] = 8;
    A[1][3] = -11;
    A[2][3] = -3;


    Matrix B = A + 3.14;

    std::cout << "Matrix A:\n" << A << std::endl;

    std::cout << "*************************" << std::endl;

    std::cout << "Matrix B = A + 3.14:\n" << B << std::endl;

    std::cout << "*************************" << std::endl;

    std::vector<double> res = A.Gauss();


    for (int i = 0; i < res.size(); i++) {
        std::cout << "x" << i + 1 << " = " << res[i] << std::endl;
    }


    Matrix C(3, 3);
    C[0][0] = 4; C[0][1] = 1; C[0][2] = 1;
    C[1][0] = 1; C[1][1] = 6; C[1][2] = -1;
    C[2][0] = 1; C[2][1] = 2; C[2][2] = 5;

    std::vector<double> y = {9, 10, 20};

    double* result = Jacobi(C, y);
    std::cout << "*************************" << std::endl;

    std::cout << "Result for Jacobi method: ";
    for (int i = 0; i < C.getRows(); ++i) {
        std::cout << result[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "*************************" << std::endl;
    delete[] result;

//    int n;
//
//    std::cout << "\nEnter the no. of data pairs to be entered:\n";
//    std::cin >> n;
//
//    Matrix x(n, 1), y(n, 1);
//
//    std::cout << "\nEnter the x-axis values:\n";
//
//    std::cin >> x;
//
//    std::cout << "\nEnter the y-axis values:\n";
//
//    std::cin >> y;
//
//    linearRegression(n, x, y);

    return 0;
}
