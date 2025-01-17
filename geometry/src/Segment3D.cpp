#include "Segment3D.h"
#include "Vector3D.h"
#include "Point3D.h"
#include "Line3D.h"
#include "Matrix.h"
#include <climits>
#include <cmath>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

Segment3D::Segment3D()
{
	m_startPoint = Point3D();
	m_endPoint = Point3D();
}

//-----------------------------------------------------------------------------

Segment3D::Segment3D(const Point3D& startPoint, const Point3D& endPoint)
{
	m_startPoint = startPoint;
	m_endPoint = endPoint;
}

//-----------------------------------------------------------------------------

Point3D Segment3D::GetStartPoint() const
{
	return m_startPoint;
}

//-----------------------------------------------------------------------------

Point3D Segment3D::GetEndPoint() const
{
	return m_endPoint;
}

//-----------------------------------------------------------------------------

void Segment3D::SetStartPoint(const Point3D& startPoint)
{
	m_startPoint = startPoint;
}

//-----------------------------------------------------------------------------

void Segment3D::SetEndPoint(const Point3D& endPoint)
{
	m_endPoint = endPoint;
}

//-----------------------------------------------------------------------------

double Segment3D::GetLength()const
{
	return sqrt(pow((m_startPoint.GetX() - m_endPoint.GetX()), 2)
		        + pow((m_startPoint.GetY() - m_endPoint.GetY()), 2) 
	            + pow((m_startPoint.GetZ() - m_endPoint.GetZ()), 2));
}

//-----------------------------------------------------------------------------

Point3D Segment3D::GetApproximateProjectionOfPointOnSegment(const Point3D& pointToProject)const
{
	Vector3D directionVectorOfLine((m_endPoint.GetX() - m_startPoint.GetX()), 
		                            (m_endPoint.GetY() - m_startPoint.GetZ()), 
		                            (m_endPoint.GetZ() - m_startPoint.GetZ()));
	std::vector<Point3D> pointsAlongVector;
	pointsAlongVector.reserve(100);

	directionVectorOfLine.GetPointsAlongVector(m_startPoint, pointsAlongVector);

	double minDistance = std::numeric_limits<double>::max();
	Point3D approximateClosestPoint;

	for (int iPoint = 0; iPoint < pointsAlongVector.size() ; ++iPoint)
	{
		Point3D currentPoint = pointsAlongVector[iPoint];
		double distanceBetweenPoints = currentPoint.GetDistance(pointToProject);
		
		if (distanceBetweenPoints < minDistance)
		{
			minDistance = distanceBetweenPoints;
			approximateClosestPoint = currentPoint;
		}
	}
	return approximateClosestPoint;
}

//-----------------------------------------------------------------------------

Point3D Segment3D::GetPreciseProjectionOfPointOnSegmentUsingVector(const Point3D& pointToProject)const
{
	Vector3D segmentVector(m_startPoint,m_endPoint);
	segmentVector.Normalize();
	Vector3D vectorToProject(m_startPoint,pointToProject);
	double lengthOfProjection = segmentVector.DotProduct(vectorToProject);

	if (lengthOfProjection >= ZeroConstant && ((lengthOfProjection-GetLength())>=ZeroConstant))
		return m_endPoint;
	else if (lengthOfProjection < ZeroConstant)
		return m_startPoint;

	return Point3D((m_startPoint.GetX() + (lengthOfProjection * segmentVector.GetI())), (m_startPoint.GetY()
		+ lengthOfProjection * segmentVector.GetJ()), (m_startPoint.GetZ() + lengthOfProjection * segmentVector.GetK()));
}

//-----------------------------------------------------------------------------

void Segment3D::GetPointsOnSegment(std::vector<Point3D>& pointsOnSegment)const
{
	Vector3D segmentAlongVector(m_startPoint,m_endPoint);
	segmentAlongVector.GetPointsAlongVector(m_startPoint, pointsOnSegment);
}

//-----------------------------------------------------------------------------

bool Segment3D::GetInstersectionOfTwoSegments(const Segment3D& otherSegment, Point3D& pointOfInterSection) const
{
	if (IsParallel(otherSegment) || IsSkew(otherSegment))
		return false;

	Vector3D directionOfFirstSegment(m_startPoint, m_endPoint);
	Vector3D directionOfSecondSegment(otherSegment.m_startPoint, otherSegment.m_endPoint);

	std::vector<std::vector<double>> coefficients(2, std::vector<double>(2, 0));
	coefficients[0].push_back(directionOfFirstSegment.GetI());
	coefficients[0].push_back(-1 * directionOfSecondSegment.GetI());

	coefficients[1].push_back(directionOfFirstSegment.GetJ());
	coefficients[1].push_back(-1 * directionOfSecondSegment.GetJ());

	std::vector<double> rightMatrix{ (otherSegment.m_startPoint.GetX() - m_startPoint.GetX()), (otherSegment.m_startPoint.GetY() - m_startPoint.GetY())};

	const double* valuesOfParameters = Matrix::SolveAXB(coefficients, rightMatrix);

	if (((fabs(valuesOfParameters[0]) >= ZeroConstant) && (fabs(valuesOfParameters[0]) <= 1))
		&& ((fabs(valuesOfParameters[1]) >= ZeroConstant) && (fabs(valuesOfParameters[1]) <= 1)))
	{
		double u = valuesOfParameters[0];
		pointOfInterSection.SetX(m_startPoint.GetX() + (u * directionOfFirstSegment.GetI()));
		pointOfInterSection.SetY(m_startPoint.GetY() + (u * directionOfFirstSegment.GetJ()));
		pointOfInterSection.SetZ(m_startPoint.GetZ() + (u * directionOfFirstSegment.GetK()));
		return true;
	}

	return false;
}

//-----------------------------------------------------------------------------

bool Segment3D::IsSkew(const Segment3D& otherSegment) const
{
	Line3D firstLine(m_startPoint,Vector3D(m_startPoint,m_endPoint));
	Line3D secondLine(otherSegment.m_startPoint,Vector3D(otherSegment.m_startPoint, otherSegment.m_endPoint));
	
	return firstLine.IsSkew(secondLine);
}

//-----------------------------------------------------------------------------

bool Segment3D::IsParallel(const Segment3D& otherSegment) const
{
	if (fabs(Vector3D(m_startPoint,m_endPoint).CalculateAngle(Vector3D(otherSegment.m_startPoint,otherSegment.m_endPoint))) <= ZeroConstant)
		return true;
	
	return false;
}

//-----------------------------------------------------------------------------

Segment3D::~Segment3D()
{

}

//-----------------------------------------------------------------------------


