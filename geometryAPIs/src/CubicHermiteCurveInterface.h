#pragma once

#include <vector>

class Point3D;
class Vector3D;

class CubicHermiteCurveInterface
{

public:

	//Function to get points along cubic hermite curve
	virtual void GetPointsAlongCubitHermiteCurve(std::vector<Point3D>& pointAlongCurve)const = 0;

	//Function to get approx point of projection on cubic Hermite curve
	virtual Point3D GetApproxProjectionPointOnCubicHermiteCurve(const Point3D& pointToProject)const = 0;

	//Default destructor
	virtual ~CubicHermiteCurveInterface() = default;

};
