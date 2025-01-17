#pragma once
#include <vector>

//Forward declarations
class Point3D;

class PolynomialInterface
{
public:

    //Function to get points along polycurve
    virtual void GetPointsAlongPolycurve(std::vector<Point3D>& pointsOnPolyCurve)const = 0;

    //Function to get approx projection point on polynomial curve
    virtual Point3D GetApproxProjectionPointOnPolynomialCurve(const Point3D& pointToProject)const = 0;

    //Default constructor
    virtual ~PolynomialInterface() = default;
};
