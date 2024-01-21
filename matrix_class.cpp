#include <iostream>
#include <vector>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows;
    int cols;

public:

    Matrix(int rows, int cols) : rows(rows), cols(cols) {
        data = std::vector<std::vector<double>>(rows, std::vector<double>(cols, 0.0));
    }

    Matrix& operator=(const Matrix& matrix2) {
        if (this != &matrix2) {
            rows = matrix2.rows;
            cols = matrix2.cols;
            data = matrix2.data;
        }
        return *this;
    }

    int getRows() const {
        return rows;
    }

    int getCols() const {
        return cols;
    }

    bool operator==(const Matrix& matrix2) const {
        return (rows == matrix2.rows && cols == matrix2.cols && data == matrix2.data);
    }

    Matrix operator+(const Matrix& matrix2) const {
        if (rows != matrix2.rows || cols != matrix2.cols) {
            throw std::invalid_argument("The dimensions of the matrices do not match");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = data[i][j] + matrix2.data[i][j];
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& matrix2) const {
        if (rows != matrix2.rows || cols != matrix2.cols) {
            throw std::invalid_argument("The dimensions of the matrices do not match");
        }

        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = data[i][j] - matrix2.data[i][j];
            }
        }
        return result;
    }

    Matrix operator*(double scalar) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = data[i][j] * scalar;
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& matrix2) const {
        if (cols != matrix2.rows) {
            throw std::invalid_argument("The number of columns in the first matrix must be equal to the number of rows in the second matrix.");
        }

        Matrix result(rows, matrix2.cols);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < matrix2.cols; j++) {
                for (int k = 0; k < cols; k++) {
                    result[i][j] += data[i][k] * matrix2.data[k][j];
                }
            }
        }

        return result;
    }

    std::vector<double> operator*(const std::vector<double>& vec) const {
        if (cols != vec.size()) {
            throw std::invalid_argument("The number of columns in the matrix must be equal to the size of the vector.");
        }

        std::vector<double> result(rows, 0.0);

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += data[i][j] * vec[j];
            }
        }

        return result;
    }

    Matrix operator~() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[j][i] = data[i][j];
            }
        }
        return result;
    }

    friend std::ostream& operator<<(std::ostream& out, const Matrix& matrix) {
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                out << matrix.data[i][j] << ' ';
            }
            out << std::endl;
        }
        return out;
    }

    friend std::istream& operator>>(std::istream& in, Matrix& matrix) {
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                in >> matrix.data[i][j];
            }
        }
        return in;
    }

    std::vector<double>& operator[](int index) {
        return data[index];
    }

    friend Matrix operator+(const Matrix& matrix, double scalar) {
        Matrix result(matrix.rows, matrix.cols);
        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                result[i][j] = matrix.data[i][j] + scalar;
            }
        }
        return result;
    }

    std::vector<double> Gauss();


};
