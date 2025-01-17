#include "Segment2D.h"
#include<cmath>


//-----------------------------------------------------------------------------

Segment2D::Segment2D()
{
	m_startPoint = Point2D();
	m_endPoint = Point2D();
}

//-----------------------------------------------------------------------------

Segment2D::Segment2D(const Point2D& startPoint, const Point2D& endPoint)
{
	m_startPoint = startPoint;
	m_endPoint = endPoint;
}

//-----------------------------------------------------------------------------

Point2D Segment2D::GetStartPoint() const
{
	return m_startPoint;
}

//-----------------------------------------------------------------------------

Point2D Segment2D::GetEndPoint() const
{
	return m_endPoint;
}

//-----------------------------------------------------------------------------

void Segment2D::SetStartPoint(const Point2D& startPoint)
{
	m_startPoint = startPoint;
}

//-----------------------------------------------------------------------------

void Segment2D::SetEndPoint(const Point2D& endPoint)
{
	m_endPoint = endPoint;
}

//-----------------------------------------------------------------------------

double Segment2D::GetLength()const
{
	return sqrt(pow((m_startPoint.GetX() - m_endPoint.GetX()),2)
					 +pow((m_startPoint.GetY() - m_endPoint.GetY()),2));
}

//-----------------------------------------------------------------------------

double Segment2D::GetSlope()const
{
	return atan((m_endPoint.GetX() - m_startPoint.GetX())/(m_endPoint.GetY() - m_startPoint.GetY()));
}

//-----------------------------------------------------------------------------

Segment2D::~Segment2D()
{
	
}

//-----------------------------------------------------------------------------


