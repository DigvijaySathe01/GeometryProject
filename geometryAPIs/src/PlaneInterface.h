#pragma once

class Line3D;
class Point3D;
class Vector3D;
class Plane;

class PlaneInterface
{

public:

	//Function to get normal
	virtual Vector3D GetNormal()const = 0;

	//Function to get point on plane
	virtual Point3D GetPointOnPlane()const = 0;

	//Function to get angle between the planes
	virtual double GetAngleBetweenPlanes(const Plane& otherPlane)const = 0;

	//Function to get distance between plane and point
	virtual double GetDistanceOfPointFromPlane(const Point3D& point)const = 0;

	//Function to get angle between plane and line
	virtual double GetAngleBetweenPlaneAndLine(const Line3D& line)const = 0;

	//Function to get projection of point on plane
	virtual Point3D GetProjectionOfPointOnPlane(const Point3D& pointToProject)const = 0;

	//Function to check whether a point is on plane or not
	virtual bool IsPointOnPlane(const Point3D& pointToCheck)const = 0;

	//Function to get intersection of line with plane
	virtual Point3D GetIntersectionOfLineWithPlane(const Line3D& line)const = 0;

	//Function to check whether two planes are parallel or not
	virtual bool IsParallel(const Plane& otherPlane)const = 0;

	//Function to get InterSectionOfTwoPlanes;
	virtual Line3D GetInterSectionOfTwoPlanes(const Plane& otherPlane)const = 0;

	//Destructor
	virtual ~PlaneInterface() = default;

};

