#include "Vector2D.h"
#include "Point2D.h"
#include<cmath>
#include<cassert>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

Vector2D::Vector2D()
{
	m_i = 0.0;
	m_j = 0.0;
}

//-----------------------------------------------------------------------------

Vector2D::Vector2D(double i, double j)
{
	m_i = i;
	m_j = j;
}

//-----------------------------------------------------------------------------

Vector2D::Vector2D(double magnitude, const Vector2D& direction)
{
	m_i = direction.m_i * magnitude;
	m_j = direction.m_j * magnitude;
}

//-----------------------------------------------------------------------------

Vector2D::Vector2D(const Point2D& startPoint, const Point2D& endPoint)
{
	m_i = (endPoint.GetX() - startPoint.GetX());
	m_j = (endPoint.GetY() - startPoint.GetY());
}

//-----------------------------------------------------------------------------

double Vector2D::GetI()const
{
	return m_i;
}

//-----------------------------------------------------------------------------

double Vector2D::GetJ()const
{
	return m_j;
}

//-----------------------------------------------------------------------------

void Vector2D::GetIJ(double& i, double& j)const
{
	i = m_i;
	j = m_j;
}

//-----------------------------------------------------------------------------

void Vector2D::SetI(double i)
{
	m_i = i;
}

//-----------------------------------------------------------------------------

void Vector2D::SetJ(double j)
{
	m_j = j;
}

//-----------------------------------------------------------------------------

void Vector2D::SetIJ(double i, double j)
{
	m_i = i;
	m_j = j;
}

//-----------------------------------------------------------------------------

Vector2D Vector2D::operator+(const Vector2D& secondVector)const
{
	return Vector2D((m_i + secondVector.m_i), (m_j + secondVector.m_j));
}

//-----------------------------------------------------------------------------

Vector2D Vector2D::operator-(const Vector2D& secondVector)const
{
	return Vector2D((m_i - secondVector.m_i), (m_j - secondVector.m_j));
}

//-----------------------------------------------------------------------------

Vector2D Vector2D::operator*(double scalarValue) const
{
	return Vector2D(scalarValue * m_i, scalarValue * m_j);
}

//-----------------------------------------------------------------------------

double Vector2D::DotProduct(const Vector2D& secondVector)const
{
	return ((m_i * secondVector.m_i) + (m_j * secondVector.m_j));
}

//-----------------------------------------------------------------------------

double Vector2D::GetMagnitude()const
{
	return sqrt((m_i * m_i) + (m_j * m_j));
}

//-----------------------------------------------------------------------------

double Vector2D::CalculateAngle(const Vector2D& secondVector)const
{
	assert(IsValid());
	assert(secondVector.IsValid());

	if (!IsValid() || !secondVector.IsValid())
		return 0;

	return acos(DotProduct(secondVector) / ((GetMagnitude()) * (secondVector.GetMagnitude())));
}

//-----------------------------------------------------------------------------

void Vector2D::Normalize()
{
	double magnitude = GetMagnitude();
	assert(IsValid());

	m_i = m_i / magnitude;
	m_j = m_j / magnitude;
}

//-----------------------------------------------------------------------------

bool Vector2D::IsValid()const
{
	if (GetMagnitude() > ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

void Vector2D::GetPointsAlongVector(const Point2D& startPoint, std::vector<Point2D>& pointsAlongVector)const
{
	for (double iParam = 0.01; iParam < 1; iParam += 0.01)
	{
		double x = startPoint.GetX() + m_i * iParam;
		double y = startPoint.GetY() + m_j * iParam;
		pointsAlongVector.emplace_back(x,y);
	}
}

//-----------------------------------------------------------------------------

Point2D Vector2D::GetPointAlongVector(const Point2D& point, double param)const
{
	return Point2D((point.GetX() + param * m_i), (point.GetY() + param * m_j));
}

//-----------------------------------------------------------------------------

Vector2D::~Vector2D()
{

}

//-----------------------------------------------------------------------------

