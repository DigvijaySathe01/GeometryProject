#pragma once

#pragma once
#include <vector>

class Matrix;

class MatrixInterface
{

public:

	//Function to get minor of the element
	virtual double GetMinor(int row, int column)const = 0;

	//Function to set elements of matrix
	virtual void SetElement(int row, int column, double value) = 0;

	//Function to get elements of matrix in single vector
	virtual const double* GetElementsOfMatrix()const = 0;

	//Function to get single element of the matrix
	virtual double GetMatrixElement(int row, int column)const = 0;

	//Function to get number of rows of matrix
	virtual int GetNumRows()const = 0;

	//Function to get number of columns of matrix
	virtual int GetNumColumns()const = 0;

	//Function to get determinant of matrix
	virtual double GetDeterminant()const = 0;

	//Function to get cofactor of the element
	virtual double GetCofactor(int row, int column)const = 0;

	//Function to get transpose of the matrix
	virtual void Transpose(Matrix& transpose)const = 0;

	//Function to inverse of the matrix
	virtual void Inverse(Matrix& inverse)const = 0;

	//Overloaded '/' function for division of matrix by scalar
	virtual void operator/(double scalar) = 0;

	//Overloaded '*' function for division of matrix by scalar
	virtual void operator*(double scalar) = 0;

	//Overloaded '+' function for addition of two matrices
	virtual Matrix operator+(const Matrix& otherMatrix)const = 0;

	//Overloaded '-' function for subtraction of two matrices
	virtual Matrix operator-(const Matrix& otherMatrix)const = 0;

	//Function to check whether the matrix is singular or not
	virtual bool IsSingular()const = 0;

	//Function to get multiplication of two matrix using formula
	virtual Matrix GetMatrixMultiplication(const Matrix& otherMatrix)const = 0;

	//Default destructor
	virtual ~MatrixInterface() = default;
};

