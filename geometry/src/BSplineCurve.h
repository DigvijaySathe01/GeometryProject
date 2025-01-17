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
	double BSplineBasisFunction(int i, int degree, double param)const override;

	//Function to get points along BSpline curve
	void GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor)const override;

	//Function to get approx point of projection on BSpline curve
	Point3D GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProjection)const override;

	//Default destructor
	~BSplineCurve();
};

