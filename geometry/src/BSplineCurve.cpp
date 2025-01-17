#include "BSplineCurve.h"
#include "Point3D.h"

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

BSplineCurve::BSplineCurve(const std::vector<Point3D>& controlPoints, int degree)
{
	m_numberOfSpans = (int)controlPoints.size() + degree;
	m_degree = degree;
	m_controlPoints = controlPoints;
	m_knotVector.resize(m_numberOfSpans + 1, 0);

	//creating uniform knot vector considering for clamped BSpline curve
	int count = 0;
	for (int iKnot = 0; iKnot <= m_numberOfSpans; ++iKnot)
	{
		if (iKnot <= m_degree)
		{
			m_knotVector[iKnot] = count;
			if (iKnot == m_degree)
				count++;

			continue;
		}

		if (iKnot + m_degree == m_numberOfSpans)
		{
			while (iKnot <= m_numberOfSpans)
			{
				m_knotVector[iKnot] = count;
				iKnot++;
			}
			break;
		}
		
			m_knotVector[iKnot] = count++;

	}
}

//-----------------------------------------------------------------------------

double BSplineCurve::BSplineBasisFunction(int i, int degree, double param) const
{
	if ((degree == 0) && ((param >= m_knotVector[i]) && (param < m_knotVector[i + 1])))
		return 1;
	else if (degree == 0)
		return 0;

	double firstComponent = 0;
	double firstComponent1 = (param - m_knotVector[i]);
	double firstComponent2 = (m_knotVector[i + degree] - m_knotVector[i]);
	if (!(fabs(firstComponent2 - ZeroConstant) <= ZeroConstant))
		firstComponent = firstComponent1 / firstComponent2;

	double secondComponent = 0;
	double secondComponent1 = (m_knotVector[i + degree + 1] - param); 
	double secondComponent2 = (m_knotVector[i + degree + 1] - m_knotVector[i + 1]);
	if (!(fabs(secondComponent2 - ZeroConstant) <= ZeroConstant))
		secondComponent = secondComponent1 / secondComponent2;

	return (firstComponent * BSplineBasisFunction(i, degree - 1, param)) + (secondComponent * BSplineBasisFunction(i + 1, degree - 1, param));
}

//-----------------------------------------------------------------------------

void BSplineCurve::GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor) const
{
	int startParam = m_knotVector[m_degree];
	int endParam = m_knotVector[m_numberOfSpans - m_degree];

	double paramIncrement = (double)(endParam - startParam) / 1000;
	for (double iParam = startParam; iParam <= endParam; iParam += paramIncrement)
	{
		double x_coordinate = 0;
		double y_coordinate = 0;
		double z_coordinate = 0;

		for (int iControlPoint = 0; iControlPoint < m_controlPoints.size(); ++iControlPoint)
		{
			double basisFunctionValue = BSplineBasisFunction(iControlPoint, m_degree, iParam);
			x_coordinate += basisFunctionValue * (m_controlPoints[iControlPoint].GetX());
			y_coordinate += basisFunctionValue * (m_controlPoints[iControlPoint].GetY());
			z_coordinate += basisFunctionValue * (m_controlPoints[iControlPoint].GetZ());
		}
		pointsAlongVecor.emplace_back(x_coordinate, y_coordinate, z_coordinate);
	}
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProject) const
{
	std::vector<Point3D> pointsAlongCurve;
	GetPointsALongBSplineCurve(pointsAlongCurve);

	double distance = std::numeric_limits<double>::max();
	Point3D projectionPoint;

	for (int iPoint = 0; iPoint < pointsAlongCurve.size(); ++iPoint)
	{
		double tempDistance = pointsAlongCurve[iPoint].GetDistance(pointToProject);
		if (tempDistance < distance)
		{
			distance = tempDistance;
			projectionPoint = pointsAlongCurve[iPoint];
		}
	}

	return projectionPoint;
}

//-----------------------------------------------------------------------------


BSplineCurve::~BSplineCurve()
{
}

//-----------------------------------------------------------------------------
