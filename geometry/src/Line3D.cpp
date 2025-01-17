#include "Line3D.h"
#include "Vector3D.h"
#include "Plane.h"
#include "Matrix.h"
#include <cassert>

//-----------------------------------------------------------------------------
const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

Line3D::Line3D(const Point3D& pointOnLine, const Vector3D& direction)
{
	m_pointOnline = pointOnLine;
	m_direction   = direction;
}

//-----------------------------------------------------------------------------

Line3D::Line3D(const Point3D& firstPoint, const Point3D& secondPoint)
{
	m_pointOnline = firstPoint;
	m_direction.SetI(firstPoint.GetX()- secondPoint.GetX());
	m_direction.SetJ(firstPoint.GetY()- secondPoint.GetY());
	m_direction.SetK(firstPoint.GetX()- secondPoint.GetZ());
}

//-----------------------------------------------------------------------------

Point3D Line3D::GetPointOnLine(double distance) const
{
	if(fabs(distance-ZeroConstant) <= ZeroConstant)
	return m_pointOnline;

	Vector3D lineDirection(m_direction);
	return lineDirection.GetPointAlongVector(m_pointOnline, 10);
}

//-----------------------------------------------------------------------------

Vector3D Line3D::GetDirectionVector() const
{
	return m_direction;
}

//-----------------------------------------------------------------------------

double Line3D::GetDistanceOfpointFromLine(const Point3D& point) const
{
	Point3D projectionOfPoint = GetProjectionOfPointOnLine(point);
	return projectionOfPoint.GetDistance(m_pointOnline);
}

//-----------------------------------------------------------------------------

Point3D Line3D::GetProjectionOfPointOnLine(const Point3D& pointToProject) const
{
	Vector3D otherVector(m_pointOnline, pointToProject);
	Vector3D tempDirection(m_direction);
	tempDirection.Normalize();
	double distance = tempDirection.DotProduct(otherVector);
	return tempDirection.GetPointAlongVector(m_pointOnline,distance);
}

//-----------------------------------------------------------------------------

bool Line3D::IsSkew(const Line3D& otherLine) const
{
	Point3D otherPointOnFirstline = GetPointOnLine(10);
	Point3D otherPointOnSecondline = otherLine.GetPointOnLine(10);

	Plane tempPlane(m_pointOnline, otherPointOnFirstline,otherLine.GetPointOnLine());

	if (tempPlane.IsPointOnPlane(otherPointOnSecondline))
		return true;

	return false;
}

//-----------------------------------------------------------------------------

double Line3D::GetDistanceBetweenSkewLines(const Line3D& otherLine) const
{
	Point3D otherPointOnFirstLine = GetPointOnLine(10);
	Point3D otherPointOnSecondtLine = otherLine.GetPointOnLine(10);

	Vector3D tempVector(otherPointOnFirstLine, otherPointOnSecondtLine);
	Vector3D CrossProductOfDirection = m_direction.CrossProduct(otherLine.GetDirectionVector());
	CrossProductOfDirection.Normalize();
	return fabs(tempVector.DotProduct(CrossProductOfDirection));
}

//-----------------------------------------------------------------------------

double Line3D::GetAngleBetweenLines(const Line3D& otherLine) const
{
	assert(!IsSkew(otherLine));
	
	return m_direction.CalculateAngle(otherLine.m_direction);
}

//-----------------------------------------------------------------------------

bool Line3D::IsParallel(const Line3D& otherLine) const
{
	if (fabs(GetAngleBetweenLines(otherLine) - ZeroConstant) <= ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

double Line3D::GetDistanceBetweenParallelLines(const Line3D& otherLine) const
{
	assert(IsParallel(otherLine));

	Vector3D tempDirection(m_direction);
	tempDirection.Normalize();
	Vector3D vectorFromTwoPontsOnTwoLines(m_pointOnline, otherLine.m_pointOnline);

	return fabs(vectorFromTwoPontsOnTwoLines.CrossProduct(tempDirection).GetMagnitude());
}

//-----------------------------------------------------------------------------

Point3D Line3D::GetItersectionOfTwoLines(const Line3D& otherLine) const
{
	assert((!IsSkew(otherLine)) && (!IsParallel(otherLine)) && "The two lines are either skew or parallel");

	Vector3D directionOfFirstLine(m_direction);
	directionOfFirstLine.Normalize();

	Vector3D directionOfSecondLine(otherLine.m_direction);
	directionOfSecondLine.Normalize();

	std::vector<std::vector<double>> coefficients(2,std::vector<double>(2,0));
	coefficients[0].push_back(directionOfFirstLine.GetI());
	coefficients[0].push_back(-1*directionOfSecondLine.GetI());

	coefficients[1].push_back(directionOfFirstLine.GetJ());
	coefficients[1].push_back(-1 * directionOfSecondLine.GetJ());

	std::vector<double> rightMatrix{(otherLine.m_pointOnline.GetX() - m_pointOnline.GetX()), (otherLine.m_pointOnline.GetY() - m_pointOnline.GetY())};

	const double* valuesOfParameters = Matrix::SolveAXB(coefficients, rightMatrix);

	assert(fabs((valuesOfParameters[0] * directionOfFirstLine.GetK()) - (valuesOfParameters[1] * directionOfSecondLine.GetK())) <= ZeroConstant);

	double u = valuesOfParameters[0];
	return Point3D(m_pointOnline.GetX() + (u * directionOfFirstLine.GetI()),
		            m_pointOnline.GetY() + (u * directionOfFirstLine.GetJ()), 
		            m_pointOnline.GetZ() + (u * directionOfFirstLine.GetK()));
}

//-----------------------------------------------------------------------------

Line3D::~Line3D()
{

}

//-----------------------------------------------------------------------------








