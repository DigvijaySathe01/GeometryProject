#pragma once
#include "Vector3DInterface.h"

//Forward declarations
class Point3D;

class Vector3D : public Vector3DInterface
{
private:
	double m_i;  //x component of vector
	double m_j;  //y component of vector
	double m_k;  //z component of vector
public:

	//Default constructor
	Vector3D();

	//Parameterized constructor with x,y and z component
	Vector3D(double i, double j, double k);

	//Parameterized constructor with start point and end point
	Vector3D(const Point3D& startPoint, const Point3D& endPoint);

	//Parameterized constructor with start point and end point
	Vector3D(double magnitude, const Vector3D& direction);

	//Parameterized constructor with point for position vector
	Vector3D(const Point3D& point);

	//Function to get x component of vector
	double GetI()const override;

	//Function to get y component of vector
	double GetJ()const override;

	//Function to get z component of vector
	double GetK()const override;

	//Function to get x, y, z components of vector
	void GetIJK(double& i, double& j, double& k)const override;

	//Function to set x component of vector
	void SetI(double i) override;

	//Function to get y component of vector
	void SetJ(double j) override;

	//Function to get z component of vector
	void SetK(double k) override;

	//Function to get x, y, z components of vector
	void SetIJK(double i, double j, double k) override;

	//Overloaded '+' operator for addition two vectors
	Vector3D operator+(const Vector3D& secondVector)const override;

	//Overloaded '-' operator for substraction of two vectors
	Vector3D operator-(const Vector3D& secondVector)const override;

	//Overloaded '*' operator for scalar Multiplication
	Vector3D operator*(double scalarValue)const override;

	//Fuction to calculage dot product of two vectors
	double DotProduct(const Vector3D& secondVector)const override;

	//Fuction to calculage cross product of two vectors
	Vector3D CrossProduct(const Vector3D& secondVector)const override;

	//Function to calculate angle between two vectors
	double CalculateAngle(const Vector3D& secondVector)const override;

	//Function to calculate magnitude of vector
	double GetMagnitude()const override;

	//Function to normalize the vector
	void Normalize() override;

	//Function to check validity of the vector
	bool IsValid()const override;

	//Function to get 100 points along vector
	void GetPointsAlongVector(const Point3D& startPoint, std::vector<Point3D>& pointsAlongVector)const override;

	//Function to get point along vector from given point
	Point3D GetPointAlongVector(const Point3D& point, double param) override;

	//Function to check whether two vectors are parallel or not
	bool isParallel(const Vector3D& otherVector)const override;

	//Function to check whether two vectors are perpendicular or not
	bool isPerpendicular(const Vector3D& otherVector)const override;

	//Destructor
	~Vector3D();
};

