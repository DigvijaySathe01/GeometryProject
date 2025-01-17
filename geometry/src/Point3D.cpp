#include "Point3D.h"
#include "Vector3D.h"
#include <cmath>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

Point3D::Point3D()
{
	m_x = 0.0;
	m_y = 0.0;
	m_z = 0.0;
}

//-----------------------------------------------------------------------------

Point3D::Point3D(double x, double y, double z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

//-----------------------------------------------------------------------------

double Point3D::GetX()const
{
	return m_x;
}

//-----------------------------------------------------------------------------

double Point3D::GetY()const
{
	return m_y;
}

//-----------------------------------------------------------------------------

double Point3D::GetZ()const
{
	return m_z;
}

//-----------------------------------------------------------------------------

void Point3D::GetXYZ(double& x, double& y, double& z)const
{
	x = m_x;
	y = m_y;
	z = m_z;
}

//-----------------------------------------------------------------------------

void Point3D::SetX(double x)
{
	m_x = x;
}

//-----------------------------------------------------------------------------

void Point3D::SetY(double y)
{
	m_y = y;
}

//-----------------------------------------------------------------------------

void Point3D::SetZ(double z)
{
	m_z = z;
}

//-----------------------------------------------------------------------------

void Point3D::SetXYZ(double x, double y, double z)
{
	m_x = x;
	m_y = y;
	m_z = z;
}

//-----------------------------------------------------------------------------

double Point3D::GetDistance(Point3D secondPoint)const
{
	return std::sqrt(std::pow((m_x - secondPoint.m_x), 2)
					+ std::pow((m_y - secondPoint.m_y), 2)
	                + std::pow((m_z - secondPoint.m_z), 2));
}

//-----------------------------------------------------------------------------

Point3D Point3D::GetPointAlongVector(const Vector3D& vector, double distance)
{
	return Point3D((m_x + distance*vector.GetI()),(m_y + distance * vector.GetJ()),(m_z + distance * vector.GetK()));
}

//-----------------------------------------------------------------------------

Point3D Point3D::operator*(const double scalar) const
{
	return Point3D(m_x * scalar, m_y * scalar, m_z * scalar);
}

//-----------------------------------------------------------------------------

Point3D Point3D::operator+(const Point3D& otherPoint) const
{
	return Point3D(m_x + otherPoint.m_x, m_y + otherPoint.m_y, m_z + otherPoint.m_z);
}

//-----------------------------------------------------------------------------

Point3D Point3D::operator-(const Point3D& otherPoint) const
{
	return Point3D(m_x - otherPoint.m_x, m_y - otherPoint.m_y, m_z - otherPoint.m_z);
}

//-----------------------------------------------------------------------------

bool Point3D::operator==(const Point3D& otherPoint) const
{
	return (fabs(m_x - otherPoint.m_x) <= ZeroConstant) && 
		   (fabs(m_y - otherPoint.m_y) <= ZeroConstant) && 
		   (fabs(m_z - otherPoint.m_z) <= ZeroConstant);
}

//-----------------------------------------------------------------------------

Point3D::~Point3D()
{
}

//-----------------------------------------------------------------------------
