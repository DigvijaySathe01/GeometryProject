#pragma once

class Point3D;
class Vector3D;

class Line3D;

class Line3DInterface
{

public:

	//Function to get point on line
	virtual Point3D GetPointOnLine(double distance = 0)const = 0;

	//Function to get direction vector of line
	virtual Vector3D GetDirectionVector()const = 0;

	//Function to get distance of point from line
	virtual double GetDistanceOfpointFromLine(const Point3D& point)const = 0;

	//Function to check whether two lines are skew are not 
	virtual bool IsSkew(const Line3D& otherLine)const = 0;

	//Function to get distance between two skew lines
	virtual double GetDistanceBetweenSkewLines(const Line3D& otherLine)const = 0;

	//Function to get angle between two lines
	virtual double GetAngleBetweenLines(const Line3D& otherLine)const = 0;

	//Function to check whether two lines are parallel or not
	virtual bool IsParallel(const Line3D& otherLine)const = 0;

	//Function to get distance between two parallel lines
	virtual double GetDistanceBetweenParallelLines(const Line3D& otherLine)const = 0;

	//Function to get projection of point on line
	virtual Point3D GetProjectionOfPointOnLine(const Point3D& pointToProject)const = 0;

	//Function to get intersection of two lines
	virtual Point3D GetItersectionOfTwoLines(const Line3D& otherLine)const = 0;

	//Destructor
	virtual ~Line3DInterface() = default;
};

