#include "BSplineCurve.h"
#include "Point3D.h"
#include <cassert>
#include <stdexcept>
#include <climits>
//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

BSplineCurve::BSplineCurve(const std::vector<Point3D>& controlPoints, int degree)
{
	assert(!controlPoints.empty() && degree > 0 && controlPoints.size() > degree && "Invalid control points or degree");

	m_numberOfSpans = static_cast<int>(controlPoints.size()) + degree;
	m_degree = degree;
	m_controlPoints = controlPoints;
	m_knotVector.resize(m_numberOfSpans + 1, 0);

	//creating uniform knot vector considering for clamped BSpline curve
	for (int iKnot = 0; iKnot <= m_numberOfSpans; ++iKnot) {
		if (iKnot <= m_degree) {
			m_knotVector[iKnot] = 0;
		}
		else if (iKnot >= m_numberOfSpans - m_degree) {
			m_knotVector[iKnot] = m_numberOfSpans - m_degree;
		}
		else {
			m_knotVector[iKnot] = iKnot - m_degree;
		}
	}
}

//-----------------------------------------------------------------------------

int BSplineCurve::GetDegree() const
{
	return m_degree;
}

//-----------------------------------------------------------------------------

const std::vector<Point3D>& BSplineCurve::GetControlPoints() const
{
	return m_controlPoints;
}

//-----------------------------------------------------------------------------

const std::vector<double>& BSplineCurve::GetKnotVector() const
{
	return m_knotVector;
}

//-----------------------------------------------------------------------------

double BSplineCurve::GetKnotAtIndex(const size_t index) const
{
	if (index >= m_knotVector.size())
		throw std::out_of_range("knot vector index out of bound!!!");

	return m_knotVector[index];
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::GetControlPointAtIndex(const size_t index) const
{
	if (index >= m_controlPoints.size())
		throw std::out_of_range("Control point vector index out of bound!!!");


	return m_controlPoints[index];
}

//-----------------------------------------------------------------------------

double BSplineCurve::BSplineBasisFunction(const int i, const int degree, const double param) const
{
	if (degree == 0)
		return ((param >= m_knotVector[i]) && (param < m_knotVector[i + 1]) ? 1.0 : 0.0);

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

Point3D BSplineCurve::GetPointOnCurve(const double param) const
{
	assert(param >= m_knotVector[m_degree] && param <= m_knotVector[m_numberOfSpans - m_degree] && "Parameter value out of bound!!!");

	double x_coordinate = 0;
	double y_coordinate = 0;
	double z_coordinate = 0;
	int controlPointsSize = static_cast<int>(m_controlPoints.size());
	for (int iControlPoint = 0; iControlPoint < controlPointsSize; ++iControlPoint)
	{
		double basisFunctionValue = BSplineBasisFunction(iControlPoint, m_degree, param);
		x_coordinate += basisFunctionValue * (m_controlPoints[iControlPoint].GetX());
		y_coordinate += basisFunctionValue * (m_controlPoints[iControlPoint].GetY());
		z_coordinate += basisFunctionValue * (m_controlPoints[iControlPoint].GetZ());
	}

	return Point3D(x_coordinate, y_coordinate, z_coordinate);
}

//-----------------------------------------------------------------------------

void BSplineCurve::GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor, const int numPoints) const
{
	assert(numPoints > 0 && "Invalid number of points!!!");

	double startParam = m_knotVector[m_degree];
	double endParam = m_knotVector[m_numberOfSpans - m_degree];

	double paramIncrement = static_cast<double>(endParam - startParam) / numPoints;
	for (double iParam = startParam; iParam <= endParam; iParam += paramIncrement)
	{
		pointsAlongVecor.emplace_back(GetPointOnCurve(iParam));
	}
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProject) const
{
	std::vector<Point3D> pointsAlongCurve;
	GetPointsALongBSplineCurve(pointsAlongCurve, 1000);

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

double BSplineCurve::BSplineBasisFunctionDerivative(const int i, const int degree, 
													const double param, const int order) const
{
	assert(param >= m_knotVector[m_degree] && param <= m_knotVector[m_numberOfSpans - m_degree] && "Parameter value out of bound!!!");

	assert(order >= 0 && order <= degree && "Invalid order of derivative");

	if (order == 0)
		return BSplineBasisFunction(i, degree, param);

	double leftTerm = 0;
	double leftNumerator = degree;
	double leftDemominator = m_knotVector[i + degree] - m_knotVector[i];
	if (fabs(leftDemominator) > ZeroConstant)
		leftTerm = (leftNumerator / leftDemominator) * BSplineBasisFunctionDerivative(i, degree-1, param, order-1);

	double rightTerm = 0;
	double rightNumerator = degree;
	double rightDemominator = m_knotVector[i + degree + 1] - m_knotVector[i + 1];
	if (fabs(rightDemominator) > ZeroConstant)
		rightTerm = (rightNumerator / rightDemominator) * BSplineBasisFunctionDerivative(i+1, degree - 1, param, order - 1);


	return leftTerm - rightTerm;
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::GetPointOnBSplineDerivative(const double param, const int order) const
{
	assert(param >= m_knotVector[m_degree] && param <= m_knotVector[m_numberOfSpans - m_degree] && "Parameter value out of bound!!!");

	assert(order >= 0 && order <= m_degree && "Invalid order of derivative");

	double xCoOrdinate = 0.0;
	double yCoOrdinate = 0.0;
	double zCoOrdinate = 0.0;
	for (int iCPoint = 0; iCPoint < m_controlPoints.size(); ++iCPoint)
	{
		double basisFunctionValue = BSplineBasisFunctionDerivative(iCPoint, m_degree, param, order);
		
		xCoOrdinate += (m_controlPoints[iCPoint].GetX() * basisFunctionValue);
		yCoOrdinate += (m_controlPoints[iCPoint].GetY() * basisFunctionValue);
		zCoOrdinate += (m_controlPoints[iCPoint].GetZ() * basisFunctionValue);
	}

	return Point3D(xCoOrdinate, yCoOrdinate, zCoOrdinate);
}

//-----------------------------------------------------------------------------

void BSplineCurve::GetControlPointsOfBSplineDerivative(const int order, std::vector<Point3D>& derivativeControlPoints) const
{
	assert(order >= 0 && order <= m_degree && "Invalid order of derivative");

	derivativeControlPoints = m_controlPoints;
	int derivativeDegree = m_degree;
	std::vector<Point3D> tempControlPoints;

	for (int iOrder = 1; iOrder <= order; ++iOrder)
	{
		tempControlPoints.clear();
		for (int iCPoint=0; iCPoint < derivativeControlPoints.size()-1; ++iCPoint)
		{
			double numerator = derivativeDegree;
			double denominator = m_knotVector[iCPoint + derivativeDegree + 1] - m_knotVector[iCPoint + 1];
			if (fabs(denominator) < ZeroConstant) {
				throw std::runtime_error("Invalid knot interval!!!");
			}

			Point3D newContolPoint = (derivativeControlPoints[iCPoint + 1] - derivativeControlPoints[iCPoint]) * (numerator / denominator);
			tempControlPoints.push_back(newContolPoint);
		}
		--derivativeDegree;
		derivativeControlPoints = tempControlPoints;
	}
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::ProjectPoint(const Point3D& pointToProject, const double guessParam, 
								   const double endParam, const int maxIterations) const
{
	return ProjectPointUsingNewtonRaphsonMethod(pointToProject, guessParam, endParam, maxIterations);
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::ProjectPointUsingBruteForceMethod(const Point3D& pointToProject, const double startParam, 
													    const double endParam, const int numSamples) const
{
	double minDistance = std::numeric_limits<double>::max();
	Point3D closestPoint;

	for (int iSample = 0; iSample < numSamples; ++iSample)
	{
		double currentParam = startParam + iSample * (endParam - startParam) / (numSamples - 1);
		Point3D currentPoint = GetPointOnCurve(currentParam);
		double distance = currentPoint.GetSquareDistance(pointToProject);
		if (distance < minDistance)
		{
			minDistance = distance;
			closestPoint = currentPoint;
		}
	}
	return closestPoint;
}

//-----------------------------------------------------------------------------

Point3D BSplineCurve::ProjectPointUsingNewtonRaphsonMethod(const Point3D& pointToProject, const double guessParam, 
														   const double endParam, const int maxIterations) const
{
	assert(guessParam >= m_knotVector[m_degree] && guessParam <= m_knotVector[m_numberOfSpans - m_degree] && "Parameter value out of bound!!!");

	double t = guessParam;
	for (int iIter = 0; iIter < maxIterations; ++iIter)
	{
		Point3D Ct = GetPointOnCurve(t);
		Point3D dCt = GetPointOnBSplineDerivative(t, 1);
		Point3D ddCt = GetPointOnBSplineDerivative(t, 2);

		Point3D diff = Point3D(Ct.GetX()-pointToProject.GetX(), Ct.GetY() - pointToProject.GetY(), Ct.GetZ() - pointToProject.GetZ());
		double firstDerivative = 2 * (diff.GetX()*dCt.GetX() + diff.GetY() * dCt.GetY() + diff.GetZ() * dCt.GetZ());
		double secondDerivative  = 2 * (dCt.GetX() * dCt.GetX() + dCt.GetY() * dCt.GetY() + dCt.GetZ() * dCt.GetZ() 
										+ diff.GetX() * ddCt.GetX() + diff.GetY() * ddCt.GetY() + diff.GetZ() * ddCt.GetZ());

		if (fabs(secondDerivative) <= ZeroConstant)
			return ProjectPointUsingBruteForceMethod(pointToProject, guessParam, endParam, 100);

		double tNew = t - firstDerivative / secondDerivative;

		tNew = std::max(guessParam,std::min(tNew,endParam));

		if (fabs(tNew - t) <= ZeroConstant)
		{
			t = tNew;
			break;
		}

		t = tNew;
	}

	return GetPointOnCurve(t);
}

//-----------------------------------------------------------------------------

