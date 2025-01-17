#include "Point2D.h"
#include "cmath"

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

 Point2D::Point2D( )
{
	m_x = 0.0;
	m_y = 0.0;
}

//-----------------------------------------------------------------------------

 Point2D::Point2D(double x, double y)
 {
	 m_x = x;
	 m_y = y;
 }

 //-----------------------------------------------------------------------------

double Point2D::GetX()const
{
	return m_x;
}

//-----------------------------------------------------------------------------


double Point2D::GetY()const
{
	return m_y;
}

//-----------------------------------------------------------------------------

void Point2D::GetXY(double& x, double& y)const
{
	x = m_x;
	y = m_y;
}

//-----------------------------------------------------------------------------

void Point2D::SetX(double x)
{
	m_x = x;
}

//-----------------------------------------------------------------------------

void Point2D::SetY(double y)
{
	m_y = y;
}

//-----------------------------------------------------------------------------

void Point2D::SetXY(double x, double y)
{
	m_x = x;
	m_y = y;
}

//-----------------------------------------------------------------------------

double Point2D::GetDistance(Point2D secondPoint)const
{
	return std::sqrt(std::pow((m_x-secondPoint.m_x),2)
					 + std::pow((m_y - secondPoint.m_y), 2));
}

//-----------------------------------------------------------------------------

Point2D Point2D::operator*(const double scalar) const
{
	return Point2D(m_x * scalar, m_y * scalar);
}

//-----------------------------------------------------------------------------

Point2D Point2D::operator+(const Point2D& otherPoint) const
{
	return Point2D(m_x + otherPoint.m_x, m_y + otherPoint.m_y);
}

//-----------------------------------------------------------------------------

Point2D Point2D::operator-(const Point2D& otherPoint) const
{
	return Point2D(m_x - otherPoint.m_x, m_y - otherPoint.m_y);
}

//-----------------------------------------------------------------------------

bool Point2D::operator==(const Point2D& otherPoint) const
{
	return (fabs(m_x - otherPoint.m_x) <= ZeroConstant) && (fabs(m_y - otherPoint.m_y) <= ZeroConstant);
}

//-----------------------------------------------------------------------------

Point2D::~Point2D( )
{

}

//-----------------------------------------------------------------------------