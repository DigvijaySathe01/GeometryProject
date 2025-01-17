#include "Polynomial.h"
#include "Matrix.h"
#include "Point3D.h"
#include <climits>

//-----------------------------------------------------------------------------

Polynomial::Polynomial(const std::vector<Point3D>& dataPointsFromUser)
{
	m_pointsGivenByUser = dataPointsFromUser;
	m_degree = (int)dataPointsFromUser.size() - 1;
	std::vector<double> coefficientVector((m_degree + 1) * (m_degree + 1), 0);

	std::vector<double> rightVectorForX(m_degree + 1, 0);
	std::vector<double> rightVectorForY(m_degree + 1, 0);
	std::vector<double> rightVectorForZ(m_degree + 1, 0);

	for (int iRow = 0; iRow <= m_degree; ++iRow)
	{
		for (int iColumn = 0; iColumn <= m_degree; ++iColumn)
		{
			coefficientVector[(iRow * (m_degree + 1)) + iColumn] = pow(iRow, iColumn);
		}
		rightVectorForX[iRow] = dataPointsFromUser[iRow].GetX();
		rightVectorForY[iRow] = dataPointsFromUser[iRow].GetY();
		rightVectorForZ[iRow] = dataPointsFromUser[iRow].GetZ();
	}

	Matrix coefficientMatrix(m_degree + 1, m_degree + 1,coefficientVector);
	Matrix inversCoefficientMatrix(coefficientMatrix.GetNumRows(), coefficientMatrix.GetNumColumns());
	coefficientMatrix.Inverse(inversCoefficientMatrix);

	Matrix righMatrixForX(m_degree + 1, 1, rightVectorForX);
	Matrix righMatrixForY(m_degree + 1, 1, rightVectorForY);
	Matrix righMatrixForZ(m_degree + 1, 1, rightVectorForZ);

	Matrix coefficientsForX = inversCoefficientMatrix.GetMatrixMultiplication(righMatrixForX);
	Matrix coefficientsForY = inversCoefficientMatrix.GetMatrixMultiplication(righMatrixForY);
	Matrix coefficientsForZ = inversCoefficientMatrix.GetMatrixMultiplication(righMatrixForZ);

	size_t size = coefficientsForX.GetNumRows() * coefficientsForX.GetNumRows();
	m_xEquationCoefficients.resize(size);
	m_yEquationCoefficients.resize(size);
	m_zEquationCoefficients.resize(size);

	for (size_t iSize = 0; iSize < size; ++iSize)
	{
		m_xEquationCoefficients[iSize] = coefficientsForX.GetMatrixElement((int)iSize, 1);
		m_yEquationCoefficients[iSize] = coefficientsForY.GetMatrixElement((int)iSize, 1);
		m_zEquationCoefficients[iSize] = coefficientsForZ.GetMatrixElement((int)iSize, 1);
	}
		
}

//-----------------------------------------------------------------------------

void Polynomial::GetPointsAlongPolycurve(std::vector<Point3D>& pointsOnPolyCurve) const
{
	double incrementFactor =0.0005;          //To get 50 points on the curve;
	
	for (double iParam = 0; iParam <= m_degree; iParam+= incrementFactor)
	{
		double xCoordinate = 0;
		double yCoordinate = 0;
		double zCoordinate = 0;

		for (int iDegree = 0; iDegree <= m_degree; ++iDegree)
		{
			xCoordinate += (m_xEquationCoefficients[iDegree] * pow(iParam, iDegree));
			yCoordinate += (m_yEquationCoefficients[iDegree] * pow(iParam, iDegree));
			zCoordinate += (m_zEquationCoefficients[iDegree] * pow(iParam, iDegree));
		}
		pointsOnPolyCurve.emplace_back(xCoordinate, yCoordinate, zCoordinate);
	}
}

//-----------------------------------------------------------------------------

Point3D Polynomial::GetApproxProjectionPointOnPolynomialCurve(const Point3D& pointToProject) const
{
	std::vector<Point3D> pointsAlongCurve;
	GetPointsAlongPolycurve(pointsAlongCurve);

	double distance = std::numeric_limits<double>::max();
	Point3D projectionPoint;

	for (int iPoint = 0; iPoint < pointsAlongCurve.size(); ++iPoint)
	{
		double tempDistance = pointsAlongCurve[iPoint].GetDistance(pointToProject);
		if (tempDistance < distance)
		{
			distance = tempDistance;
			projectionPoint = pointsAlongCurve[iPoint];
		}
	}

	return projectionPoint;
}

//-----------------------------------------------------------------------------

Polynomial::~Polynomial()
{

}

//-----------------------------------------------------------------------------

