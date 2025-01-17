#pragma once
#include "PolynomialInterface.h"

//Forward declarations
class Point3D;

class Polynomial : public PolynomialInterface
{
private:
    int m_degree;                                  // Degree of the polynomial
    std::vector<double> m_xEquationCoefficients;   // Coefficients of the parametric polynomial equation for the x coordinate
    std::vector<double> m_yEquationCoefficients;   // Coefficients of the parametric polynomial equation for the y coordinate
    std::vector<double> m_zEquationCoefficients;   // Coefficients of the parametric polynomial equation for the z coordinate
    std::vector<Point3D> m_pointsGivenByUser;      // Vector to store points given by user

public:

    //Parameterized constructor with data points from user
    Polynomial(const std::vector<Point3D>& dataPointsFromUser);

    //Function to get points along polycurve
    void GetPointsAlongPolycurve(std::vector<Point3D>& pointsOnPolyCurve)const override;

    //Function to get approx projection point on polynomial curve
    Point3D GetApproxProjectionPointOnPolynomialCurve(const Point3D& pointToProject)const override;

    //Default constructor
    ~Polynomial();
};
