#pragma once
#include "Point3D.h"
#include "Vector3D.h"

#include "PlaneInterface.h"

class Line3D;

class Plane : public PlaneInterface
{
private:
	Point3D m_pointOnPlane;      //Point on plane
	Vector3D m_normal;			  //Normal of plane

public:

	//Parameterized constructor with three points
	Plane(const Point3D& firstPoint, const Point3D& secondPoint, const Point3D& thirdPoint);

	//Parameterized constructor with point on plane and normal vector
	Plane(const Point3D& pointOnPlane, const Vector3D& normal);

	//Parameterized constructor with two vectors on plane and point on plane
	Plane(const Vector3D& firstVectorOnPlane, const Vector3D& secondVectorOnPlane, Point3D& pointOnPlane);

	//Parameterized constructor with plane and distance from plane
	Plane(const Plane& otherPlane, double distance);

	//Function to get normal
	Vector3D GetNormal()const override;

	//Function to get point on plane
	Point3D GetPointOnPlane()const override;

	//Function to get angle between the planes
	double GetAngleBetweenPlanes(const Plane& otherPlane)const override;

	//Function to get distance between plane and point
	double GetDistanceOfPointFromPlane(const Point3D& point)const override;

	//Function to get angle between plane and line
	double GetAngleBetweenPlaneAndLine(const Line3D& line)const override;

	//Function to get projection of point on plane
	Point3D GetProjectionOfPointOnPlane(const Point3D& pointToProject)const override;

	//Function to check whether a point is on plane or not
	bool IsPointOnPlane(const Point3D& pointToCheck)const override;

	//Function to get intersection of line with plane
	Point3D GetIntersectionOfLineWithPlane(const Line3D& line)const override;

	//Function to check whether two planes are parallel or not
	bool IsParallel(const Plane& otherPlane)const override;

	//Function to get InterSectionOfTwoPlanes;
	Line3D GetInterSectionOfTwoPlanes(const Plane& otherPlane)const override;

	//Destructor
	~Plane();

};

