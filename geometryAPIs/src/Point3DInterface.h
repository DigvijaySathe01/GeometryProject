#pragma once

//Forward declarations
class Vector3D;
class Point3D;

class Point3DInterface 
{
public:

	//Function to get x co-ordinate of the point 
	virtual double GetX()const = 0;

	//Function to get y co-ordinate of the point
	virtual double GetY()const = 0;

	//Function to get z co-ordinate of the point
	virtual double GetZ()const = 0;

	//Function to get x, y and z co-ordinate of the point
	virtual void  GetXYZ(double& x, double& y, double& z)const = 0;

	//Function to set x co-ordinate of the point 
	virtual void SetX(double x) = 0;

	//Function to set y co-ordinate of the point 
	virtual void SetY(double y) = 0;

	//Function to set z co-ordinate of the point
	virtual void SetZ(double z) = 0;

	//Function to set x,y and z co-ordinate of the point 
	virtual void SetXYZ(double x, double y, double z) = 0;

	//This function calculate the distance between two points
	virtual double GetDistance(Point3D secondPoint)const = 0;

	//Function to get point along vector
	virtual Point3D GetPointAlongVector(const Vector3D& vector, double distance) = 0;

	//Overloaded * operator to scale a Point3D by a scalar value
	virtual Point3D operator*(const double scalar)const = 0;

	//Overloaded + operator to add two points
	virtual Point3D operator+(const Point3D& otherPoint)const = 0;

	//Overloaded - operator to add two points
	virtual Point3D operator-(const Point3D& otherPoint)const = 0;

	//Overloaded == operator to check equality of two points
	virtual bool operator==(const Point3D& otherPoint)const = 0;

	//Destructor
	virtual ~Point3DInterface() = default;
};