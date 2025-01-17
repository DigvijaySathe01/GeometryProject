#pragma once

//Standard header
#include <vector>
#include <memory>

//Forward declarations
class Point3D;

class BezierCurveInterface
{

public:

    //Function to get points along bazier curve
    virtual void GetPoinstAlongBazierCurve(const int numPoints, std::vector<Point3D>& pointsAlongCurve)const = 0;

    //Function to get approx point of projection on Bezier curve
    virtual Point3D GetApproxPointOfProjectionOnBezierCurve(const Point3D& pointToProject)const = 0;

    //Function to insert a control point at given index
    virtual void AddControlPoint(const size_t index, const Point3D& controlPoint) = 0;

    //Function to remove control point
    virtual void RemoveControlPoint(const size_t index) = 0;

    //Function to update the control point at given index
    virtual void UpdadeControlPoint(const size_t index, const Point3D& controlPoint) = 0;

    //Function to get control point at given index
    virtual void GetControlPoint(const size_t index, Point3D& controlPoint)const = 0;

    //Function to get all control points of Bezier curve
    virtual void GetControlPoints(std::vector<Point3D>& controlPoints)const = 0;

    //Function to set control points of the curve
    virtual void SetControlPoints(const std::vector<Point3D>& controlPoints) = 0;

    //Function to get degree of the curve
    virtual int GetDegree()const = 0;

    //Function to get number of control points
    virtual size_t GetNumControlPoints()const = 0;

    //Function to get point along curve for given parametric value
    virtual void GetPointOnCurve(const double param, Point3D& pointOnCurve)const = 0;

    //Function to get point along curve using De Casteljau's algorithm
    virtual void GetPointOnCurveUsingDeCasteljau(const double param, Point3D& point)const = 0;

    //Function to split curve at given parametric value
    virtual void SplitCurve(const double param, std::vector<std::unique_ptr<BezierCurveInterface>>& splittedCurves)const = 0;

    //function to compute n th order derivative
    virtual void GetNthOrderDerivative(const int order, std::vector<Point3D>& derivativeControlPoints)const = 0;

    //Function to check C2 parametric continuity of the curve
    virtual bool IsC2Continuous(const double param)const = 0;

    //Default destructor
    virtual ~BezierCurveInterface() = default;
};