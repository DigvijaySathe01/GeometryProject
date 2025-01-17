#pragma once

//Forward declarations
class Point2D;

class Segment2DInterface
{
public:

	//This function returns the start point of the line
	virtual Point2D GetStartPoint()const = 0;

	//This function returns the end point of the line
	virtual Point2D GetEndPoint()const = 0;

	//Function to set start point of the line
	virtual void SetStartPoint(const Point2D& startPoint) = 0;

	//Function to set end point of the line
	virtual void SetEndPoint(const Point2D& endPoint) = 0;

	//This function returns the length of the line
	virtual double GetLength()const = 0;

	//Function to calculate slope of the line
	virtual double GetSlope()const = 0;

	//Destructor
	virtual ~Segment2DInterface() = default;
};

