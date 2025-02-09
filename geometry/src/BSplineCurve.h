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
	std::vector<double> m_knotVector;             //knot vector 

public:

	//Parameterized constructor with control points and degree of the curve
	BSplineCurve(const std::vector<Point3D>& controlPoints, int degree);

	//Function to get the degree of BSpline curve
	int GetDegree()const override;

	//Function to get control points of BSpline curve
	const std::vector<Point3D>&  GetControlPoints()const override;

	//Function to get knot vector of BSpline curve
	const std::vector<double>& GetKnotVector()const override;

	//Function to get knot at give index
	double GetKnotAtIndex(const size_t index)const override;

	//Function to get control point at given index
	Point3D GetControlPointAtIndex(const size_t index)const override;

	//Function for basis function for BSpline curve
	double BSplineBasisFunction(const int i, const int degree, const double param)const override;

	//Function to get point on BSpline curve for given parameter
	virtual Point3D GetPointOnCurve(const double param)const;

	//Function to get points along BSpline curve
	void GetPointsALongBSplineCurve(std::vector<Point3D>& pointsAlongVecor, const int numPoints)const override;

	//Function to get approx point of projection on BSpline curve
	Point3D GetApproxProjectionPointOnBSplineCurve(const Point3D& pointToProjection)const override;

	//Funtion to get n the order derivative of BSpline basis function
	double BSplineBasisFunctionDerivative(const int i, const int degree, const double param, const int order)const override;

	//Function to get point on n th order BSpline derivative for given value of paramter
	Point3D GetPointOnBSplineDerivative(const double param, const int order)const override;

	//Function to get control points of n th order BSpline derivative
	void GetControlPointsOfBSplineDerivative(const int order, std::vector<Point3D>& derivativeControlPoints)const override;

	//Function to project point on BSpline curve
	virtual Point3D ProjectPoint(const Point3D& pointToProject, const double guessParam, const double endParam, const int maxIterations = 50)const override;

	//Function to get projection of point on BSpline using Brute Force Method
	virtual Point3D ProjectPointUsingBruteForceMethod(const Point3D& pointToProject, const double startParam, const double endParam, const int numSamples)const override;

	//Function to get projection of point on BSpline using Newton-Raphson's Method
	virtual Point3D ProjectPointUsingNewtonRaphsonMethod(const Point3D& pointToProject, const double guessParam, const double endParam, const int maxInterations)const override;

	//Function to insert a knot in BSpline curve
	void InsertKnot(const double newKnot) override;

	//Function to split the BSpline curve at given knot
	void SplitCurveAtKnot(const double knot, std::vector<BSplineCurve>& splitCurve)const override;

	//Function to check the validity of the BSpline curve
	bool IsValid()const override;

	//Default destructor
	~BSplineCurve() = default;
};

