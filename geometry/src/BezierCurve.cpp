#include "BezierCurve.h"
#include "Point3D.h"
#include "Matrix.h"
#include "Vector3D.h"

//Standard includes
#include <cassert>

//-----------------------------------------------------------------------------

const double ZeroConstant = 1e-6;

//-----------------------------------------------------------------------------

BezierCurve::BezierCurve(const std::vector<Point3D>& cotrolPoints) :m_degree((int)cotrolPoints.size() - 1),
m_controlPoints(cotrolPoints)
{
}

//-----------------------------------------------------------------------------

int BezierCurve::factorial(int number) const
{
	if (number == 0 || number == 1)
		return 1;

	return number * factorial(number - 1);
}

//-----------------------------------------------------------------------------

Point3D BezierCurve::deCasteljauAlgorithm(const double param, const int i, const int j) const
{
	if (i == 0)
		return m_controlPoints[j-1];

	return (deCasteljauAlgorithm(param, i - 1, j) * (1 - param)) +
		   (deCasteljauAlgorithm(param, i - 1, j + 1) * param);
}

//-----------------------------------------------------------------------------

void BezierCurve::GetPoinstAlongBazierCurve(const int numPoints, std::vector<Point3D>& pointsAlongCurve) const
{
	double param = (double)1 / numPoints;        //To get 1000 points along the curve
	for (double iParam = 0; iParam <= 1; iParam = iParam + param)
	{
		pointsAlongCurve.push_back(deCasteljauAlgorithm(iParam, m_degree, 1));
	}
}

//-----------------------------------------------------------------------------

Point3D BezierCurve::GetApproxPointOfProjectionOnBezierCurve(const Point3D& pointToProject) const
{
	std::vector<Point3D> pointsAlongCurve;
	GetPoinstAlongBazierCurve(1000, pointsAlongCurve);

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

void BezierCurve::AddControlPoint(const size_t index, const Point3D& controlPoint)
{
	assert(index <= m_controlPoints.size() && "Invalid Index!");

	m_controlPoints.insert(m_controlPoints.begin() + index, controlPoint);
	m_degree += 1;
}

//-----------------------------------------------------------------------------

void BezierCurve::RemoveControlPoint(const size_t index)
{
	assert(index < m_controlPoints.size() && "Invalid Index!");

	m_controlPoints.erase(m_controlPoints.begin() + index);
	m_degree -= 1;
}

//-----------------------------------------------------------------------------

void BezierCurve::UpdadeControlPoint(const size_t index, const Point3D& controlPoint)
{
	assert(index < m_controlPoints.size() && "Invalid Index!");

	m_controlPoints[index] = controlPoint;
}

//-----------------------------------------------------------------------------

void BezierCurve::GetControlPoint(const size_t index, Point3D& controlPoint) const
{
	assert(index < m_controlPoints.size() && "Invalid Index!");

	controlPoint = m_controlPoints[index];
}

//-----------------------------------------------------------------------------

void BezierCurve::GetControlPoints(std::vector<Point3D>& controlPoints) const
{
	controlPoints = m_controlPoints;
}

//-----------------------------------------------------------------------------

void BezierCurve::SetControlPoints(const std::vector<Point3D>& controlPoints)
{
	assert(controlPoints.size() >= 2 && "Invalid input for curve!");

	m_controlPoints = controlPoints;
	m_degree = (int)m_controlPoints.size() - 1;
}

//-----------------------------------------------------------------------------

int BezierCurve::GetDegree() const
{
	return m_degree;
}

//-----------------------------------------------------------------------------

size_t BezierCurve::GetNumControlPoints() const
{
	return m_controlPoints.size();
}

//-----------------------------------------------------------------------------

void BezierCurve::GetPointOnCurve(const double param, Point3D& pointOnCurve) const
{
	double xCoordinate = 0.0;
	double yCoordinate = 0.0;
	double zCoordinate = 0.0;

	for (int iDegree = 0; iDegree <= m_degree; ++iDegree)
	{
		double basisFunction = (factorial(m_degree) / (factorial(m_degree - iDegree) * factorial(iDegree)))
			* pow(param, iDegree) * pow((1 - param), m_degree - iDegree);

		xCoordinate += basisFunction * m_controlPoints[iDegree].GetX();
		yCoordinate += basisFunction * m_controlPoints[iDegree].GetY();
		zCoordinate += basisFunction * m_controlPoints[iDegree].GetZ();
	}
	
	pointOnCurve.SetXYZ(xCoordinate, yCoordinate, zCoordinate);
}

//-----------------------------------------------------------------------------

void BezierCurve::GetPointOnCurveUsingDeCasteljau(const double param, Point3D& point) const
{
	point =  deCasteljauAlgorithm(param, m_degree, 1);
}

//-----------------------------------------------------------------------------

void BezierCurve::SplitCurve(const double param, std::vector<std::unique_ptr<BezierCurveInterface>>& splittedCurves) const
{
	Point3D commonControlPoint = deCasteljauAlgorithm(param, m_degree, 1);

	std::vector<Point3D> leftCurveControlPoints{ commonControlPoint };
	std::vector<Point3D> rightCurveControlPoints{ commonControlPoint };

	int iLeft = m_degree-1;
	int jLeft = 1;

	int iRight = m_degree-1;
	int jRight = 2;

	for (int index = 0; index < m_degree; ++index)
	{
		leftCurveControlPoints.emplace_back(deCasteljauAlgorithm(param, iLeft--, jLeft));
		rightCurveControlPoints.emplace_back(deCasteljauAlgorithm(param, iRight--, jRight++));
	}

	reverse(leftCurveControlPoints.begin(), leftCurveControlPoints.end());

	splittedCurves.emplace_back(std::make_unique<BezierCurve>(leftCurveControlPoints));
	splittedCurves.emplace_back(std::make_unique<BezierCurve>(rightCurveControlPoints));
}

//-----------------------------------------------------------------------------

void BezierCurve::GetNthOrderDerivative(const int order, std::vector<Point3D>& derivativeControlPoints) const
{
	if (order > m_degree || order < 0)
		return;

	derivativeControlPoints = m_controlPoints;

	ComputeDerivativeControlPoints(order, derivativeControlPoints);
}

//-----------------------------------------------------------------------------

bool BezierCurve::IsC2Continuous(const double param) const
{
	std::vector<std::unique_ptr<BezierCurveInterface>> splitCurves;
	SplitCurve(param, splitCurves);

	if (splitCurves.empty() || splitCurves.size())
	{
		assert("Bezier curve didn't get split");
		return false;
	}

	std::vector<Point3D> leftCurveCP;
	splitCurves[0]->GetControlPoints(leftCurveCP);

	std::vector<Point3D> rightCurveCP;
	splitCurves[0]->GetControlPoints(rightCurveCP);

	if (leftCurveCP.empty() || rightCurveCP.empty())
		return false;

	size_t leftCurveDegree = leftCurveCP.size() - 1;
	size_t rightCurveDegree = rightCurveCP.size() - 1;

	if (leftCurveDegree != m_degree || rightCurveDegree != m_degree)
		return false;

	//To check C0 continuity
	if (!(leftCurveCP[leftCurveDegree] == rightCurveCP[0]))
		return false;

	//To check C1 continuity 
	Vector3D leftCurveTangent(leftCurveCP[leftCurveDegree] - leftCurveCP[leftCurveDegree - 1]);
	Vector3D rightCurveTangent(rightCurveCP[1] - rightCurveCP[0]);

	if (!(fabs(leftCurveTangent.GetMagnitude() - rightCurveTangent.GetMagnitude()) <= ZeroConstant && leftCurveTangent.isParallel(rightCurveTangent)))
		return false;

	//To check C2 continuity 
	std::vector<Point3D> leftCurveDerivativeCP;
	std::vector<Point3D> rightCurveDerivativeCP;

	splitCurves[0]->GetNthOrderDerivative(2, leftCurveDerivativeCP);
	splitCurves[1]->GetNthOrderDerivative(2, rightCurveDerivativeCP);

	if (leftCurveDerivativeCP.empty() || rightCurveDerivativeCP.empty())
		return false;

	if (!(leftCurveDerivativeCP[leftCurveDerivativeCP.size() - 1] == rightCurveDerivativeCP[0]))
		return false;

	return true;
}


//-----------------------------------------------------------------------------

void BezierCurve::ComputeDerivativeControlPoints(const int order, std::vector<Point3D>& derivativeControlPoints) const
{
	if (order == 0)
		return;

	std::vector<Point3D> newControlPoints;
	size_t degree = derivativeControlPoints.size() - 1;

	for (size_t iCPoint = 0; iCPoint < degree; ++iCPoint)
	{
		Point3D derivativeControlPoint = (derivativeControlPoints[iCPoint + 1] - derivativeControlPoints[iCPoint]) * (double)degree;
		newControlPoints.push_back(derivativeControlPoint);
	}

	ComputeDerivativeControlPoints(order - 1, newControlPoints);
	derivativeControlPoints = newControlPoints;
}

//-----------------------------------------------------------------------------

BezierCurve::~BezierCurve()
{

}

//-----------------------------------------------------------------------------
