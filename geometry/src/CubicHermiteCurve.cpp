#include "CubicHermiteCurve.h"


//-----------------------------------------------------------------------------

CubicHermiteCurve::CubicHermiteCurve(const std::vector<std::pair<Point3D, Vector3D>>& pointsWithTangents)
{
	m_pointsOnCurveWithTangents = pointsWithTangents;
}

//-----------------------------------------------------------------------------

void CubicHermiteCurve::GetPointsAlongCubitHermiteCurve(std::vector<Point3D>& pointAlongCurve) const
{
	for (int iPoint = 1; iPoint < m_pointsOnCurveWithTangents.size(); ++iPoint)
	{
		double param = (double)1 / 1000;
		for (double iParam = 0; iParam <= 1; iParam += param)
		{
			double h0 = 2 * pow(iParam, 3) - 3 * pow(iParam, 2) + 1;       //h0(t) = 2t^3 - 3t^2 + 1
			double h1 = pow(iParam, 3) - 2 * pow(iParam, 2) + iParam;      //h1(t) = t^3 - 2t^2 + t
			double h2 = -2 * pow(iParam, 3) + 3 * pow(iParam, 2);          //h2(t) = -2t^3 + 3t^2
			double h3 = pow(iParam, 3) - pow(iParam, 2);                   //h3(t) = t^3 - t^2

			double xCoeff = (h0 * m_pointsOnCurveWithTangents[iPoint-1].first.GetX()) 
				            + (h1 * m_pointsOnCurveWithTangents[iPoint - 1].second.GetI())
				            + (h2 * m_pointsOnCurveWithTangents[iPoint].first.GetX()) 
				            + (h3 * m_pointsOnCurveWithTangents[iPoint].second.GetI());

			double YCoeff = (h0 * m_pointsOnCurveWithTangents[iPoint - 1].first.GetY())
							+ (h1 * m_pointsOnCurveWithTangents[iPoint - 1].second.GetJ())
							+ (h2 * m_pointsOnCurveWithTangents[iPoint].first.GetY())
							+ (h3 * m_pointsOnCurveWithTangents[iPoint].second.GetJ());

			double ZCoeff = (h0 * m_pointsOnCurveWithTangents[iPoint - 1].first.GetZ())
				            + (h1 * m_pointsOnCurveWithTangents[iPoint - 1].second.GetK())
							+ (h2 * m_pointsOnCurveWithTangents[iPoint].first.GetZ())
							+ (h3 * m_pointsOnCurveWithTangents[iPoint].second.GetK());

			pointAlongCurve.emplace_back(Point3D(xCoeff, YCoeff, ZCoeff));
		}
	}
}

//-----------------------------------------------------------------------------

Point3D CubicHermiteCurve::GetApproxProjectionPointOnCubicHermiteCurve(const Point3D& pointToProject) const
{
	std::vector<Point3D> pointsAlongCurve;
	GetPointsAlongCubitHermiteCurve(pointsAlongCurve);

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

CubicHermiteCurve::~CubicHermiteCurve()
{

}

//-----------------------------------------------------------------------------
