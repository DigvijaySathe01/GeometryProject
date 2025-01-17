#pragma once
#include "Segment3DInterface.h"
#include "Point3D.h"

class Segment3D : public Segment3DInterface
{
private:
	Point3D m_startPoint;   //start point of the line
	Point3D m_endPoint;	 //end point of the line

public:

	//Default constructor
	Segment3D();

	//Parameterized constructor
	Segment3D(const Point3D& startPoint, const Point3D& endPoint);

	//This function returns the start point of the line
	Point3D GetStartPoint()const override;

	//This function returns the end point of the line
	Point3D GetEndPoint()const override;

	//Function to set start point of the line
	void SetStartPoint(const Point3D& startPoint) override;

	//Function to set end point of the line
	void SetEndPoint(const Point3D& endPoint) override;

	//This function returns the length of the line
	double GetLength()const override;

	//Function to get approximate projection of point on line
	Point3D GetApproximateProjectionOfPointOnSegment(const Point3D& pointToProject)const override;

	//Function to get precise projection of point on line
	Point3D GetPreciseProjectionOfPointOnSegmentUsingVector(const Point3D& pointToProject)const override;

	//Function to get points along segment
	void GetPointsOnSegment(std::vector<Point3D>& pointsOnSegment)const override;

	//Function to get intersection of two planes
	bool GetInstersectionOfTwoSegments(const Segment3D& otherSegment, Point3D& pointOfInterSection)const override;

	//Function to check whether two segments are skew or not
	bool IsSkew(const Segment3D& otherSegment)const override;

	//Function to check whether two segments are parallel or not
	bool IsParallel(const Segment3D& otherSegment)const override;

	//Destructor
	~Segment3D();
};

