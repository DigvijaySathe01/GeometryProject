#include "Vector3D.h"
#include "Point3D.h"
#include<cmath>
#include<cassert>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

Vector3D::Vector3D()
{
	m_i = 0.0;
	m_j = 0.0;
	m_k = 0.0;
}

//-----------------------------------------------------------------------------

Vector3D::Vector3D(double i, double j, double k)
{
	m_i = i;
	m_j = j;
	m_k = k;
}

//-----------------------------------------------------------------------------

Vector3D::Vector3D(double magnitude, const Vector3D& direction)
{
	m_i = direction.m_i * magnitude;
	m_j = direction.m_j * magnitude;
	m_k = direction.m_k * magnitude;
}

//-----------------------------------------------------------------------------

Vector3D::Vector3D(const Point3D& startPoint, const Point3D& endPoint)
{
	m_i = (endPoint.GetX() - startPoint.GetX());
	m_j = (endPoint.GetY() - startPoint.GetY());
	m_k = (endPoint.GetZ() - startPoint.GetZ());
}

//-----------------------------------------------------------------------------

Vector3D::Vector3D(const Point3D& point)
{
	m_i = point.GetX();
	m_j = point.GetY();
	m_k = point.GetZ();
}

//-----------------------------------------------------------------------------


double Vector3D::GetI()const
{
	return m_i;
}

//-----------------------------------------------------------------------------

double Vector3D::GetJ()const
{
	return m_j;
}

//-----------------------------------------------------------------------------

double Vector3D::GetK()const
{
	return m_k;
}

//-----------------------------------------------------------------------------

void Vector3D::GetIJK(double& i, double& j, double& k)const
{
	i = m_i;
	j = m_j;
	k = m_k;
}

//-----------------------------------------------------------------------------

void Vector3D::SetI(double i)
{
	m_i = i;
}

//-----------------------------------------------------------------------------

void Vector3D::SetJ(double j)
{
	m_j = j;
}

//-----------------------------------------------------------------------------

void Vector3D::SetK(double k)
{
	m_k = k;
}

//-----------------------------------------------------------------------------

void Vector3D::SetIJK(double i, double j, double k)
{
	m_i = i;
	m_j = j;
	m_k = k;
}

//-----------------------------------------------------------------------------

Vector3D Vector3D::operator+(const Vector3D& secondVector)const
{
	return Vector3D((m_i+secondVector.m_i), (m_j + secondVector.m_j), (m_k + secondVector.m_k));
}

//-----------------------------------------------------------------------------

Vector3D Vector3D::operator-(const Vector3D& secondVector)const
{
	return Vector3D((m_i - secondVector.m_i), (m_j - secondVector.m_j), (m_k - secondVector.m_k));
}

//-----------------------------------------------------------------------------

Vector3D Vector3D::operator*(double scalarValue) const
{
	return Vector3D(scalarValue*m_i, scalarValue * m_j, scalarValue * m_k);
}

//-----------------------------------------------------------------------------

double Vector3D::DotProduct(const Vector3D& secondVector)const
{
	return ((m_i*secondVector.m_i) + (m_j * secondVector.m_j) + (m_k * secondVector.m_k));
}

//-----------------------------------------------------------------------------

Vector3D Vector3D::CrossProduct(const Vector3D& secondVector)const
{
	double xComponent = m_i*((m_j*secondVector.m_k)-(m_k*secondVector.m_j));
	double yComponent = -m_j*((m_i*secondVector.m_k)-(m_k*secondVector.m_i));
	double zComponent = m_k*((m_i*secondVector.m_j)-(m_j*secondVector.m_i));
	return Vector3D(xComponent, yComponent, zComponent);
}

//-----------------------------------------------------------------------------

double Vector3D::GetMagnitude()const
{
	return sqrt((m_i*m_i) + (m_j*m_j) + (m_k*m_k));
}

//-----------------------------------------------------------------------------

double Vector3D::CalculateAngle(const Vector3D& secondVector)const
{
	assert(IsValid());
	assert(secondVector.IsValid());

	if (!IsValid() || !secondVector.IsValid())
		return 0;

	return acos(DotProduct(secondVector)/((GetMagnitude())*(secondVector.GetMagnitude())));
}

//-----------------------------------------------------------------------------

void Vector3D::Normalize()
{
	double magnitude = GetMagnitude();
	assert(IsValid());

	m_i = m_i / magnitude;
	m_j = m_j / magnitude;
	m_k = m_k / magnitude;

}

//-----------------------------------------------------------------------------

bool Vector3D::IsValid()const
{
	if (GetMagnitude() > ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

void Vector3D::GetPointsAlongVector(const Point3D& startPoint, std::vector<Point3D>& pointsAlongVector)const
{
	for (double iParam = 0.01 ; iParam < 1 ; iParam+=0.01)
	{
		Point3D pointOnVector((startPoint.GetX() + iParam*m_i), (startPoint.GetY() + iParam * m_j), (startPoint.GetZ() + iParam * m_k));
		pointsAlongVector.push_back(pointOnVector);
	}
}

//-----------------------------------------------------------------------------

Point3D Vector3D::GetPointAlongVector(const Point3D& point, double param)
{
	Normalize();
	return Point3D((point.GetX() + param * m_i), (point.GetY() + param * m_j), (point.GetZ() + param * m_k));
}

//-----------------------------------------------------------------------------

bool Vector3D::isPerpendicular(const Vector3D& otherVector)const
{
	double pi = 3.1415926535897;
	if (fabs(pi/2 - CalculateAngle(otherVector)) <= ZeroConstant)
		return true;
	return false;
}

//-----------------------------------------------------------------------------

bool Vector3D::isParallel(const Vector3D& otherVector)const
{
	if (fabs(CalculateAngle(otherVector)) <= ZeroConstant)
		return true;

	return false;
}

//-----------------------------------------------------------------------------

Vector3D::~Vector3D()
{

}

//-----------------------------------------------------------------------------

