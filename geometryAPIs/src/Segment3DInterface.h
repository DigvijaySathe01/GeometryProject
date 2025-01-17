#pragma once
#include <vector>

//Forward declarations
class Point3D;
class Segment3D;

class Segment3DInterface
{
public:

	//This function returns the start point of the line
	virtual Point3D GetStartPoint()const = 0;

	//This function returns the end point of the line
	virtual Point3D GetEndPoint()const = 0;

	//Function to set start point of the line
	virtual void SetStartPoint(const Point3D& startPoint) = 0;

	//Function to set end point of the line
	virtual void SetEndPoint(const Point3D& endPoint) = 0;

	//This function returns the length of the line
	virtual double GetLength()const = 0;

	//Function to get approximate projection of point on line
	virtual Point3D GetApproximateProjectionOfPointOnSegment(const Point3D& pointToProject)const = 0;

	//Function to get precise projection of point on line
	virtual Point3D GetPreciseProjectionOfPointOnSegmentUsingVector(const Point3D& pointToProject)const = 0;

	//Function to get points along segment
	virtual void GetPointsOnSegment(std::vector<Point3D>& pointsOnSegment)const = 0;

	//Function to get intersection of two planes
	virtual bool GetInstersectionOfTwoSegments(const Segment3D& otherSegment, Point3D& pointOfInterSection)const = 0;

	//Function to check whether two segments are skew or not
	virtual bool IsSkew(const Segment3D& otherSegment)const = 0;

	//Function to check whether two segments are parallel or not
	virtual bool IsParallel(const Segment3D& otherSegment)const = 0;

	//Destructor
	virtual ~Segment3DInterface() = default;
};

