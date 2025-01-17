#pragma once
#include <vector>

//Forward Declaration
class Point3D;
class Vector3D;

class Vector3DInterface
{
public:

	//Function to get x component of vector
	virtual double GetI()const = 0;

	//Function to get y component of vector
	virtual double GetJ()const = 0;

	//Function to get z component of vector
	virtual double GetK()const = 0;

	//Function to get x, y, z components of vector
	virtual void GetIJK(double& i, double& j, double& k)const = 0;

	//Function to set x component of vector
	virtual void SetI(double i) = 0;

	//Function to get y component of vector
	virtual void SetJ(double j) = 0;

	//Function to get z component of vector
	virtual void SetK(double k) = 0;

	//Function to get x, y, z components of vector
	virtual void SetIJK(double i, double j, double k) = 0;

	//Overloaded '+' operator for addition two vectors
	virtual Vector3D operator+(const Vector3D& secondVector)const = 0;

	//Overloaded '-' operator for substraction of two vectors
	virtual Vector3D operator-(const Vector3D& secondVector)const = 0;

	//Overloaded '*' operator for scalar Multiplication
	virtual Vector3D operator*(double scalarValue)const = 0;

	//Fuction to calculage dot product of two vectors
	virtual double DotProduct(const Vector3D& secondVector)const = 0;

	//Fuction to calculage cross product of two vectors
	virtual Vector3D CrossProduct(const Vector3D& secondVector)const = 0;

	//Function to calculate angle between two vectors
	virtual double CalculateAngle(const Vector3D& secondVector)const = 0;

	//Function to calculate magnitude of vector
	virtual double GetMagnitude()const = 0;

	//Function to normalize the vector
	virtual void Normalize() = 0;

	//Function to check validity of the vector
	virtual bool IsValid()const = 0;

	//Function to get 100 points along vector
	virtual void GetPointsAlongVector(const Point3D& startPoint, std::vector<Point3D>& pointsAlongVector)const = 0;

	//Function to get point along vector from given point
	virtual Point3D GetPointAlongVector(const Point3D& point, double param) = 0;

	//Function to check whether two vectors are parallel or not
	virtual bool isParallel(const Vector3D& otherVector)const = 0;

	//Function to check whether two vectors are perpendicular or not
	virtual bool isPerpendicular(const Vector3D& otherVector)const = 0;

	//Destructor
	virtual ~Vector3DInterface() = default;
};

