#pragma once
#include "PolyLine3DInterface.h"

//Forward declarations
class Segment3D;
class Point3D;

class PolyLine3D : public PolyLine3DInterface
{
private:
	std::vector<Segment3D*> m_segments;

public:

	//Parameterized constructor
	PolyLine3D(const std::vector<Point3D>& pointsFromLine);

	//Function to get 100 points on nth segment of polyline
	void GetPointsAlongSegment(std::vector<Point3D>& pointsOnSegments, int index) override;

	//Function to get points on the segments of polyline
	void GetPointsAlongSegments(std::vector<Point3D>& pointsOnSegments) override;

	//Function to get precise point of intersection on polyline
	Point3D GetPrecisePointOfProjection(const Point3D& pointToProject)const override;

	//Function to get precise point of intersection on polyline
	Point3D GetApproximateOfProjection(const Point3D& pointToProject)const override;

	//Destructor
	~PolyLine3D();
};

