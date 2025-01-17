#pragma once

#include "MatrixInterface.h"
class Matrix : public MatrixInterface
{
private:
	double* m_elements;     //Pointer pointing to the dynamically allocated 1D array of matrix elements
	int m_rows;				//Number of rows of matrix
	int m_columns;			//Number of columns of matrix

public:

	//Parameterized constructor with number of rows and columns
	Matrix(int rows, int columns);

	//Parameterized constructor with number of rows and columns
	Matrix(int rows, int columns, const std::vector<double>& elements);

	//Copy constructor
	Matrix(const Matrix& otherMatrix);

	//Function to get minor of the element
	double GetMinor(int row, int column)const;

	//Function to set elements of matrix
	void SetElement(int row, int column, double value);

	//Function to get elements of matrix in single vector
	const double* GetElementsOfMatrix()const;

	//Function to get single element of the matrix
	double GetMatrixElement(int row, int column)const;

	//Function to get number of rows of matrix
	int GetNumRows()const override;

	//Function to get number of columns of matrix
	int GetNumColumns()const override;

	//Function to get determinant of matrix
	double GetDeterminant()const override;

	//Function to get cofactor of the element
	double GetCofactor(int row, int column)const override;

	//Function to get transpose of the matrix
	void Transpose(Matrix& transpose)const override;

	//Function to inverse of the matrix
	void Inverse(Matrix& inverse)const override;

	//Overloaded '/' function for division of matrix by scalar
	void operator/(double scalar) override;

	//Overloaded '*' function for division of matrix by scalar
	void operator*(double scalar) override;

	//Overloaded '+' function for addition of two matrices
	Matrix operator+(const Matrix& otherMatrix)const override;

	//Overloaded '-' function for subtraction of two matrices
	Matrix operator-(const Matrix& otherMatrix)const override;

	//Function to check whether the matrix is singular or not
	bool IsSingular()const override;

	//Function to get multiplication of two matrix using formula
	Matrix GetMatrixMultiplication(const Matrix& otherMatrix)const override;

	// Solves the linear equation system Ax = B and returns the solution vector x.
	static const double* SolveAXB(const std::vector<std::vector<double>>& coefficients, const std::vector<double>& rightMatrix);

	//Default destructor
	~Matrix();
};

