#pragma once
#include <vector>

//Forward declarations
class Point2D;
class Vector2D;

class Vector2DInterface
{
public:

	//Function to get x component of vector
	virtual double GetI()const = 0;

	//Function to get y component of vector
	virtual double GetJ()const = 0;

	//Function to get x, y, z components of vector
	virtual void GetIJ(double& i, double& j)const = 0;

	//Function to set x component of vector
	virtual void SetI(double i) = 0;

	//Function to get y component of vector
	virtual void SetJ(double j) = 0;

	//Function to get x, y, z components of vector
	virtual void SetIJ(double i, double j) = 0;

	//Overloaded '+' operator for addition two vectors
	virtual Vector2D operator+(const Vector2D& secondVector)const = 0;

	//Overloaded '-' operator for substraction of two vectors
	virtual Vector2D operator-(const Vector2D& secondVector)const = 0;

	//Overloaded '*' operator for scalar Multiplication
	virtual Vector2D operator*(double scalarValue)const = 0;

	//Fuction to calculage dot product of two vectors
	virtual double DotProduct(const Vector2D& secondVector)const = 0;

	//Function to calculate angle between two vectors
	virtual double CalculateAngle(const Vector2D& secondVector)const = 0;

	//Function to calculate magnitude of vector
	virtual double GetMagnitude()const = 0;

	//Function to normalize the vector
	virtual void Normalize() = 0;

	//Function to check validity of the vector
	virtual bool IsValid()const = 0;

	//Function to get 100 points along vector
	virtual void GetPointsAlongVector(const Point2D& startPoint, std::vector<Point2D>& pointsAlongVector)const = 0;

	//Function to get point along vector from given point
	virtual Point2D GetPointAlongVector(const Point2D& point, double param)const = 0;

	//Destructor
	virtual ~Vector2DInterface() = default;
};

