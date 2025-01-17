#include "Matrix.h"
#include <cassert>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------


Matrix::Matrix(int rows, int columns)
{
	assert(!rows && !columns && "number of rows and number of columns should not be zero");
	m_rows = rows;
	m_columns = columns;
	m_elements = new double[m_rows * m_columns] {};
}

//-----------------------------------------------------------------------------

Matrix::Matrix(int rows, int columns, const std::vector<double>& elements)
{
	m_rows = rows;
	m_columns = columns;


	size_t size = elements.size();
	m_elements = new double[size] {};

	for (int index = 0; index < size; ++index)
	{
		m_elements[index] = elements[index];
	}

}

//-----------------------------------------------------------------------------

Matrix::Matrix(const Matrix& otherMatrix)
{
	m_rows = otherMatrix.m_rows;
	m_columns = otherMatrix.m_columns;

	int size = m_rows * m_columns;
	m_elements = new double[size] {};

	for (int index = 0; index < size; ++index)
	{
		m_elements[index] = otherMatrix.m_elements[index];
	}
}

//-----------------------------------------------------------------------------

double Matrix::GetMinor(int row, int column) const
{
	Matrix minor(m_rows - 1, m_columns - 1);

	int localRow = 0, localColumn = 0;
	for (int iRow = 0; iRow < m_rows; ++iRow)
	{
		localColumn = 0;
		if (iRow == row - 1)
			continue;

		for (int iColumn = 0; iColumn < m_columns; ++iColumn)
		{
			if (iColumn == column - 1)
				continue;

			minor.m_elements[(localRow * m_columns) + iColumn] = m_elements[(iRow* m_columns) + iColumn];
			localColumn++;
		}
		localRow++;
	}
	return minor.GetDeterminant();
}

//-----------------------------------------------------------------------------

void Matrix::SetElement(int row, int column, double value)
{
	assert((row <= m_rows && row >= 0) && (column <= m_columns && row >= 0) && "index of row of column vector is out of bound");

	m_elements[(row - 1) + ((row - 1) * m_columns) + (column - 1) - 1] = value;
}

//-----------------------------------------------------------------------------

const double* Matrix::GetElementsOfMatrix() const
{
	return m_elements;
}

//-----------------------------------------------------------------------------

double Matrix::GetMatrixElement(int row, int column) const
{
	assert((row <= m_rows && row >= 0) && (column <= m_columns && row >= 0) && "index of row of column vector is out of bound");

	return m_elements[(row-1) + ((row-1) * m_columns) + (column - 1) - 1];
}

//-----------------------------------------------------------------------------

int Matrix::GetNumRows() const
{
	return m_rows;
}

//-----------------------------------------------------------------------------

int Matrix::GetNumColumns() const
{
	return m_columns;
}

//-----------------------------------------------------------------------------   

double Matrix::GetDeterminant() const
{
	assert(m_rows == m_columns && "The number of rows and the number of columns are not the same.");

	if (m_rows == 1 && m_columns == 1)
		return m_elements[0];

	double determinant = 0;
	for (int iColumn = 0; iColumn < m_columns; ++iColumn)
	{
		determinant += m_elements[iColumn] * GetCofactor(1, iColumn + 1);
	}
	return determinant;
}

//-----------------------------------------------------------------------------

double Matrix::GetCofactor(int row, int column)const
{
	return pow(-1, (row + column)) * GetMinor(row, column);
}

//-----------------------------------------------------------------------------

void Matrix::Transpose(Matrix& transpose) const
{
	transpose.m_rows = m_columns;
	transpose.m_columns = m_rows;
	delete transpose.m_elements;
	transpose.m_elements = new double[m_rows * m_columns];

	for (int iRow = 0; iRow < m_rows; ++iRow)
	{
		for (int iColumn = 0; iColumn < m_columns; ++iColumn)
		{
			transpose.m_elements[iRow * m_columns + iColumn] = m_elements[iColumn * m_columns + iRow];
		}
	}
}

//-----------------------------------------------------------------------------    

void Matrix::Inverse(Matrix& inverse) const
{
	assert((m_rows == m_columns) && "Not a square matrix");

	Matrix tempMatrix(*this);
	inverse.m_rows = m_rows;
	inverse.m_columns = m_columns;
	delete inverse.m_elements;
	inverse.m_elements = new double[m_rows * m_columns];


	for (int iRow = 0; iRow < m_rows; ++iRow)
	{
		for (int iColumn = 0; iColumn < m_columns; ++iColumn)
		{
			if (iRow == iColumn)
			{
				inverse.m_elements[(iRow * m_columns) + iRow] = 1;
			}
			else
				inverse.m_elements[(iRow * m_columns) + iColumn] = 0;
		}
	}

	for (int iRow = 0; iRow < m_rows; ++iRow)
	{
		double tempValue = tempMatrix.m_elements[(iRow * m_columns) + iRow];
		for (int iColumn = 0; iColumn < m_columns; ++iColumn)
		{
			tempMatrix.m_elements[(iRow * m_columns) + iColumn] /= tempValue;
			inverse.m_elements[(iRow * m_columns) + iColumn] /= tempValue;

		}
		for (int jRow = 0; jRow < m_rows; ++jRow)
		{
			double tempValue1 = tempMatrix.m_elements[(jRow * m_columns) + iRow];
			if (iRow == jRow)
				continue;

			for (int jColumn = 0; jColumn < m_columns; ++jColumn)
			{
				inverse.m_elements[(jRow * m_columns) + jColumn] += (-1) * inverse.m_elements[(iRow * m_columns) + jColumn] * tempValue1;
				tempMatrix.m_elements[(jRow * m_columns) + jColumn] += (-1) * tempMatrix.m_elements[(iRow * m_columns) + jColumn] * tempValue1;
			}
		}
	}
}

//-----------------------------------------------------------------------------

void Matrix::operator/(double scalar)
{
	return (this)->operator*(1 / scalar);
}

//-----------------------------------------------------------------------------

void Matrix::operator*(double scalar)
{
	int size = m_rows * m_columns;
	for (int index = 0; index < size; ++index)
	{
		m_elements[index] *= scalar;
	}
}

//-----------------------------------------------------------------------------

Matrix Matrix::operator+(const Matrix& otherMatrix) const
{
	assert(m_rows == (otherMatrix.m_rows) && (m_columns == otherMatrix.m_columns) && "Either number of row or number columns is not same");

	Matrix additionMatrix(m_rows, m_columns);

	int size = m_rows * m_columns;
	for (int index = 0; index < size; ++index)
	{
		additionMatrix.m_elements[index] = otherMatrix.m_elements[index] + m_elements[index];
	}

	return additionMatrix;
}

//-----------------------------------------------------------------------------

Matrix Matrix::operator-(const Matrix& otherMatrix) const
{
	assert(m_rows == (otherMatrix.m_rows) && (m_columns == otherMatrix.m_columns) && "Either number of row or number columns is not same");

	Matrix subtractionMatrix(m_rows, m_columns);

	int size = m_rows * m_columns;
	for (int index = 0; index < size; ++index)
	{
		subtractionMatrix.m_elements[index] = otherMatrix.m_elements[index] - m_elements[index];
	}

	return subtractionMatrix;
}

//-----------------------------------------------------------------------------

bool Matrix::IsSingular() const
{
	if (fabs(GetDeterminant() - ZeroConstant) <= ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

Matrix Matrix::GetMatrixMultiplication(const Matrix& otherMatrix) const
{
	Matrix multiplication(m_rows, otherMatrix.m_columns);

	for (int iNewRow = 0; iNewRow < multiplication.m_rows; ++iNewRow)
	{
		for (int iNewColumn = 0; iNewColumn < multiplication.m_columns; ++iNewColumn)
		{

			for (int iFirstColumn = 0; iFirstColumn < m_columns; ++iFirstColumn)
			{
				multiplication.m_elements[iNewRow + (iNewRow * m_columns) + iNewColumn - 1] +=
					m_elements[(iNewRow * m_columns) + iFirstColumn] * otherMatrix.m_elements[(iFirstColumn * m_columns) + iNewColumn];

			}
		}
	}
	return multiplication;
}

//-----------------------------------------------------------------------------

const double* Matrix::SolveAXB(const std::vector<std::vector<double>>& coefficients,
	const std::vector<double>& rightMatrix)
{
	assert(coefficients.size() == coefficients[0].size() && "coefficient matrix is not a square matrix");

	Matrix coefficientMatrix((int)coefficients.size(), (int)coefficients[0].size());
	for (int iRow = 0; iRow < coefficientMatrix.m_rows; ++iRow)
	{
		for (int iColumn = 0; iColumn < coefficientMatrix.m_columns; ++iColumn)
		{
			coefficientMatrix.m_elements[(iRow * coefficientMatrix.m_columns) + iColumn] = coefficientMatrix.m_elements[(iRow * coefficientMatrix.m_columns) + iColumn];
		}
	}

	Matrix rightHandMatrix((int)rightMatrix.size(), 1);
	for (int iRow = 0; iRow < rightHandMatrix.m_rows; ++iRow)
	{
		rightHandMatrix.m_elements[(iRow * rightHandMatrix.m_columns) + 0] = rightMatrix[iRow];
	}

	Matrix inverseOfCoefficientMatrix(coefficientMatrix.m_rows, coefficientMatrix.m_columns);
	coefficientMatrix.Inverse(inverseOfCoefficientMatrix);

	Matrix solution = inverseOfCoefficientMatrix.GetMatrixMultiplication(rightHandMatrix);

	return solution.m_elements;
}

//-----------------------------------------------------------------------------

Matrix::~Matrix()
{
	delete[] m_elements;
}

//-----------------------------------------------------------------------------







