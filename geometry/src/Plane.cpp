#include "Plane.h"
#include "Line3D.h"
#include <cassert>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

Plane::Plane(const Point3D& firstPoint, const Point3D& secondPoint, const Point3D& thirdPoint)
{
	m_pointOnPlane = firstPoint;
	Vector3D firstVector(firstPoint, secondPoint);
	Vector3D secondVector(firstPoint, thirdPoint);
	m_normal = firstVector.CrossProduct(secondVector);
}

//-----------------------------------------------------------------------------

Plane::Plane(const Point3D& pointOnPlane, const Vector3D& normal)
{
	m_pointOnPlane = pointOnPlane;
	m_normal = normal;
}

//-----------------------------------------------------------------------------

Plane::Plane(const Vector3D& firstVectorOnPlane, const Vector3D& secondVectorOnPlane, Point3D& pointOnPlane)
{
	m_pointOnPlane = pointOnPlane;
	m_normal = firstVectorOnPlane.CrossProduct(secondVectorOnPlane);
}

//-----------------------------------------------------------------------------

Plane::Plane(const Plane& otherPlane, double distance)
{
	m_normal = otherPlane.m_normal;
	m_normal.GetPointAlongVector(otherPlane.m_pointOnPlane,distance);
}

//-----------------------------------------------------------------------------

Vector3D Plane::GetNormal() const
{
	return m_normal;
}

//-----------------------------------------------------------------------------

Point3D Plane::GetPointOnPlane() const
{
	return m_pointOnPlane;
}

//-----------------------------------------------------------------------------


double Plane::GetAngleBetweenPlanes(const Plane& otherPlane) const
{
	return m_normal.CalculateAngle(otherPlane.m_normal);
}

//-----------------------------------------------------------------------------

double Plane::GetDistanceOfPointFromPlane(const Point3D& point) const
{
	Vector3D positionVectorOfPointOnPlane(m_pointOnPlane);
	Vector3D positionVectorOfPointOutsidePlane(point);
	Vector3D tempVector(m_normal);
	tempVector.Normalize();
	return fabs(tempVector.DotProduct(positionVectorOfPointOnPlane) - tempVector.DotProduct(positionVectorOfPointOutsidePlane));
}

//-----------------------------------------------------------------------------

Point3D Plane::GetProjectionOfPointOnPlane(const Point3D& pointToProject)const
{
	double distanceOfPointFromPlane = GetDistanceOfPointFromPlane(pointToProject);
	Vector3D tempVector(m_normal);
	tempVector.Normalize();

	return Point3D((pointToProject.GetX()- (distanceOfPointFromPlane * tempVector.GetI())),
		            (pointToProject.GetY() - (distanceOfPointFromPlane * tempVector.GetJ())),
		            (pointToProject.GetZ() - (distanceOfPointFromPlane * tempVector.GetK())));
}

//-----------------------------------------------------------------------------

double Plane::GetAngleBetweenPlaneAndLine(const Line3D& line) const
{
	double pi = 3.1415926535897;
	return fabs((pi/2) - m_normal.CalculateAngle(line.GetDirectionVector()));
}

//-----------------------------------------------------------------------------

bool Plane::IsPointOnPlane(const Point3D& pointToCheck) const
{
	Vector3D tempVector(m_pointOnPlane, pointToCheck);
	if (tempVector.DotProduct(m_normal) <= ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

Point3D Plane::GetIntersectionOfLineWithPlane(const Line3D& line) const
{
	Vector3D lineDirection = line.GetDirectionVector();
	lineDirection.Normalize();

	Vector3D tempNormal(m_normal);
	tempNormal.Normalize();


	Vector3D positionVectorOfPointOnPlane(m_pointOnPlane);
	Point3D pointOnLine(line.GetPointOnLine());
	Vector3D positionVectorOfPointOnLine(pointOnLine);

	double param = ((tempNormal.DotProduct(positionVectorOfPointOnPlane)) - (tempNormal.DotProduct(positionVectorOfPointOnLine))) / (tempNormal.DotProduct(lineDirection));

	return Point3D(pointOnLine.GetX() + (param* lineDirection.GetI()), pointOnLine.GetY() + (param * lineDirection.GetJ()), pointOnLine.GetZ() + (param * lineDirection.GetK()));

}

//-----------------------------------------------------------------------------

bool Plane::IsParallel(const Plane& otherPlane) const
{
	if (fabs(GetAngleBetweenPlanes(otherPlane) - ZeroConstant) <= ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

Line3D Plane::GetInterSectionOfTwoPlanes(const Plane& otherPlane) const
{
	assert(!IsParallel(otherPlane) && "Intersection not possible...two planes are parallel");

	Vector3D directionOfIntersectingLine = otherPlane.m_normal.CrossProduct(m_normal);

	Vector3D tempVector = otherPlane.m_normal.CrossProduct(directionOfIntersectingLine);

	Line3D lineInterSectingWithFirstPlane(otherPlane.m_pointOnPlane,tempVector);

	Point3D pointOnLineOfIntersection = GetIntersectionOfLineWithPlane(lineInterSectingWithFirstPlane);

	return Line3D(pointOnLineOfIntersection,directionOfIntersectingLine);
}

//-----------------------------------------------------------------------------


Plane::~Plane()
{
}

//-----------------------------------------------------------------------------

