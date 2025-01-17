#pragma once

#include <vector>

class Point3D;

class BSplineCurveInterface
{
public:

	//Function for basis function for BSpline curve
	virtual double BSplineBasisFunction(int i, int degree, double param)const = 0;

	//Function to get points along BSpline curve
	virtual void GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor)const = 0;

	//Function to get approx point of projection on BSpline curve
	virtual Point3D GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProjection)const = 0;

	//Default destructor
	virtual ~BSplineCurveInterface() = default;
};

