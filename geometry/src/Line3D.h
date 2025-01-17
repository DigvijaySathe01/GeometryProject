#pragma once
#include "Point3D.h"
#include "Vector3D.h"

#include "Line3DInterface.h"

class Line3D : public Line3DInterface
{
private:
	Point3D m_pointOnline;   //Point on line
	Vector3D m_direction;	  //Direction of line

public:

	//Parameterized constructor with point on line and direction vector
	Line3D(const Point3D& pointOnLine, const Vector3D& direction);

	//Parameterized constructor with point on line and direction vector
	Line3D(const Point3D& firstPoint, const Point3D& secondPoint);

	//Function to get point on line
	Point3D GetPointOnLine(double distance = 0)const override;

	//Function to get direction vector of line
	Vector3D GetDirectionVector()const override;

	//Function to get distance of point from line
	double GetDistanceOfpointFromLine(const Point3D& point)const override;

	//Function to check whether two lines are skew are not 
	bool IsSkew(const Line3D& otherLine)const override;

	//Function to get distance between two skew lines
	double GetDistanceBetweenSkewLines(const Line3D& otherLine)const override;

	//Function to get angle between two lines
	double GetAngleBetweenLines(const Line3D& otherLine)const override;

	//Function to check whether two lines are parallel or not
	bool IsParallel(const Line3D& otherLine)const override;

	//Function to get distance between two parallel lines
	double GetDistanceBetweenParallelLines(const Line3D& otherLine)const override;

	//Function to get projection of point on line
	Point3D GetProjectionOfPointOnLine(const Point3D& pointToProject)const override;

	//Function to get intersection of two lines
	Point3D GetItersectionOfTwoLines(const Line3D& otherLine)const override;

	//Destructor
	~Line3D();
};

