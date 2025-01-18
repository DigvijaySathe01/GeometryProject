#pragma once

//Geom header
#include "BSplineCurveInterface.h"
class Point3D;

class BSplineCurve : public BSplineCurveInterface
{
private:
	int m_numberOfSpans;                       //Number of knot spans
	int m_degree;							   //Degree of bspline curve
	std::vector<Point3D> m_controlPoints;	   //Control points of the curve
	std::vector<int> m_knotVector;             //knot vector 

public:

	//Parameterized constructor with control points and degree of the curve
	BSplineCurve(const std::vector<Point3D>& controlPoints, int degree);

	//Function for basis function for BSpline curve
	double BSplineBasisFunction(const int i, const int degree, const double param)const override;

	//Function to get points along BSpline curve
	void GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor, const int numPoints)const override;

	//Function to get approx point of projection on BSpline curve
	Point3D GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProjection)const override;

	//Funtion to get n the order derivative of BSpline basis function
	double BSplineBasisFunctionDerivative(const int i, const int degree, const double param, const int order)const;

	//Function to get point on n th order BSpline derivative for given value of paramter
	Point3D GetPointOnBSplineDerivative(const double param, const int order)const;

	//Function to get control points of n th order BSpline derivative
	void GetControlPointsOfBSplineDerivative(const int order, std::vector<Point3D>& derivativeControlPoints)const;

	//Default destructor
	~BSplineCurve() = default;
};

