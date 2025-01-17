#pragma once
#include <vector>

//Forward declarations
class Segment3D;
class Point3D;

class PolyLine3DInterface
{
public:

	//Function to get 100 points on nth segment of polyline
	virtual void GetPointsAlongSegment(std::vector<Point3D>& pointsOnSegments, int index) = 0;

	//Function to get points on the segments of polyline
	virtual void GetPointsAlongSegments(std::vector<Point3D>& pointsOnSegments) = 0;

	//Function to get precise point of intersection on polyline
	virtual Point3D GetPrecisePointOfProjection(const Point3D& pointToProject)const = 0;

	//Function to get precise point of intersection on polyline
	virtual Point3D GetApproximateOfProjection(const Point3D& pointToProject)const = 0;

	//Destructor
	virtual ~PolyLine3DInterface() = default;
};

