#pragma once

#include <vector>

class Point3D;

class BSplineCurveInterface
{
public:

	//Function for basis function for BSpline curve
	virtual double BSplineBasisFunction(const int i, const int degree, const double param)const = 0;

	//Function to get point on BSpline curve for given parameter
	virtual Point3D GetPointOnCurve(const double param)const = 0;

	//Function to get points along BSpline curve
	virtual void GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor, const int numPoints)const = 0;

	//Function to get approx point of projection on BSpline curve
	virtual Point3D GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProjection)const = 0;

	//Funtion to get n the order derivative of BSpline basis function
	virtual double BSplineBasisFunctionDerivative(const int i, const int degree, const double param, const int order)const = 0;

	//Function to get point on n th order BSpline derivative for given value of paramter
	virtual Point3D GetPointOnBSplineDerivative(const double param, const int order)const = 0;

	//Function to get control points of n th order BSpline derivative
	virtual void GetControlPointsOfBSplineDerivative(const int order, std::vector<Point3D>& derivativeControlPoints)const = 0;

	//Function to project point on BSpline curve
	virtual Point3D ProjectPoint(const Point3D& pointToProject, const double guessParam, const double endParam, const int maxIterations = 50)const = 0;

	//Function to get projection of point on BSpline using Brute Force Method
	virtual Point3D ProjectPointUsingBruteForceMethod(const Point3D& pointToProject, const double startParam, const double endParam, const int numSamples)const = 0;

	//Function to get projection of point on BSpline using Newton-Raphson's Method
	virtual Point3D ProjectPointUsingNewtonRaphsonMethod(const Point3D& pointToProject, const double guessParam, const double endParam, const int maxInterations)const = 0;

	//Default destructor
	virtual ~BSplineCurveInterface() = default;
};

