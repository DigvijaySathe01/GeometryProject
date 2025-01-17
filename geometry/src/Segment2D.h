#pragma once
#include "Segment2DInterface.h"
#include"Point2D.h"

class Segment2D : public Segment2DInterface
{
private:
	Point2D m_startPoint;   //start point of the line
	Point2D m_endPoint;	 //end point of the line

public:

	//Default constructor
	Segment2D();

	//Parameterized constructor
	Segment2D(const Point2D& startPoint, const Point2D& endPoint);

	//This function returns the start point of the line
	Point2D GetStartPoint()const override;

	//This function returns the end point of the line
	Point2D GetEndPoint()const override;

	//Function to set start point of the line
	void SetStartPoint(const Point2D& startPoint) override;

	//Function to set end point of the line
	void SetEndPoint(const Point2D& endPoint) override;

	//This function returns the length of the line
	double GetLength()const override;

	//Function to calculate slope of the line
	double GetSlope()const override;

	//Destructor
	~Segment2D();
};

