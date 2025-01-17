#pragma once
#include "Point3D.h"
#include "Vector3D.h"

#include "CubicHermiteCurveInterface.h"

class CubicHermiteCurve : public CubicHermiteCurveInterface
{
private:
	std::vector<std::pair<Point3D, Vector3D>> m_pointsOnCurveWithTangents;

public:

	//Parameterized constructor with points and their respective tanget
	CubicHermiteCurve(const std::vector<std::pair<Point3D, Vector3D>>& pointsWithTangents);

	//Function to get points along cubic hermite curve
	void GetPointsAlongCubitHermiteCurve(std::vector<Point3D>& pointAlongCurve)const;

	//Function to get approx point of projection on cubic Hermite curve
	Point3D GetApproxProjectionPointOnCubicHermiteCurve(const Point3D& pointToProject)const;

	//Default destructor
	~CubicHermiteCurve();
};

