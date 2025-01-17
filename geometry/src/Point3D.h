#pragma once
#include "Point3DInterface.h"

//Forward declarations
class Vector3D;

class Point3D : public Point3DInterface
{
private:
	double m_x;    //x co-ordinate of point
	double m_y;	   //y co-ordinate of point
	double m_z;    //z co-ordinate of point

public:

	//Default constructor
	Point3D();

	//Parameterized constructor
	Point3D(double x, double y, double z);

	//Function to get x co-ordinate of the point 
	double GetX()const override;

	//Function to get y co-ordinate of the point
	double GetY()const override;

	//Function to get z co-ordinate of the point
	double GetZ()const override;

	//Function to get x, y and z co-ordinate of the point
	void  GetXYZ(double& x, double& y, double& z)const override;

	//Function to set x co-ordinate of the point 
	void SetX(double x) override;

	//Function to set y co-ordinate of the point 
	void SetY(double y) override;

	//Function to set z co-ordinate of the point
	void SetZ(double z) override;

	//Function to set x,y and z co-ordinate of the point 
	void SetXYZ(double x, double y, double z) override;

	//This function calculate the distance between two points
	double GetDistance(Point3D secondPoint)const override;

	//Function to get point along vector
	Point3D GetPointAlongVector(const Vector3D& vector, double distance) override;

	//Overloaded * operator to scale a Point3D by a scalar value
	Point3D operator*(const double scalar)const override;

	//Overloaded + operator to add two points
	virtual Point3D operator+(const Point3D& otherPoint)const override;

	//Overloaded - operator to add two points
	virtual Point3D operator-(const Point3D& otherPoint)const override;

	//Overloaded == operator to check equality of two points
	virtual bool operator==(const Point3D& otherPoint)const override;

	//Destructor
	~Point3D();
};