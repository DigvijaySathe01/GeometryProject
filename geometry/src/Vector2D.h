#pragma once
#include "Vector2DInterface.h"
#include <vector>

//Forward declarations
class Point2D;

class Vector2D : public Vector2DInterface
{
private:
	double m_i;  //x component of vector
	double m_j;  //y component of vector
	
public:

	//Default constructor
	Vector2D();

	//Parameterized constructor with x,y component
	Vector2D(double i, double j);

	//Parameterized constructor with start point and end point
	Vector2D(const Point2D& startPoint, const Point2D& endPoint);

	//Parameterized constructor with start point and end point
	Vector2D(double magnitude, const Vector2D& direction);

	//Function to get x component of vector
	double GetI()const override;

	//Function to get y component of vector
	double GetJ()const override;

	//Function to get x, y, z components of vector
	void GetIJ(double& i, double& j)const override;

	//Function to set x component of vector
	void SetI(double i) override;

	//Function to get y component of vector
	void SetJ(double j) override;

	//Function to get x, y, z components of vector
	void SetIJ(double i, double j) override;

	//Overloaded '+' operator for addition two vectors
	Vector2D operator+(const Vector2D& secondVector)const override;

	//Overloaded '-' operator for substraction of two vectors
	Vector2D operator-(const Vector2D& secondVector)const override;

	//Overloaded '*' operator for scalar Multiplication
	Vector2D operator*(double scalarValue)const override;

	//Fuction to calculage dot product of two vectors
	double DotProduct(const Vector2D& secondVector)const override;

	//Function to calculate angle between two vectors
	double CalculateAngle(const Vector2D& secondVector)const override;

	//Function to calculate magnitude of vector
	double GetMagnitude()const override;

	//Function to normalize the vector
	void Normalize() override;

	//Function to check validity of the vector
	bool IsValid()const override;

	//Function to get 100 points along vector
	void GetPointsAlongVector(const Point2D& startPoint, std::vector<Point2D>& pointsAlongVector)const override;

	//Function to get point along vector from given point
	Point2D GetPointAlongVector(const Point2D& point, double param)const override;

	//Destructor
	~Vector2D();
};

