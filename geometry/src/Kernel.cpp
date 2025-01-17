#include <iostream>
#include "FileHandling.h"
#include "Point3D.h"
#include "PolyLine3D.h"
#include "Segment3D.h"
#include "Matrix.h"
#include "Polynomial.h"
#include "BezierCurve.h"
#include "CubicHermiteCurve.h"
#include "BSplineCurve.h"
#include <string>

int main()
{
	//int row, column;
	//std::cout<<"ënter number of row : \n";
	//std::cin >> row;

	//	std::cout << "ënter number of column : \n";
	//std::cin >> column;

	//Matrix temp(row, column);
	//

	//for (int i = 0;i< row;i++)
	//{
	//	for (int j=0; j< column;j++)
	//	{
	//		double num = 0;
	//		std::cout << "enter the number\n";
	//		std::cin >> num;
	//		temp.SetElement(i+1,j+1,num);
	//    }
	//}

	///*Matrix inverse(row,column);
	//temp.Inverse(inverse);
	//Matrix transpose(row,column);
	//temp.Transpose(transpose);*/
	//Matrix second(row,column);

	//for (int i = 0; i < row; i++)
	//{
	//	for (int j = 0; j < column; j++)
	//	{
	//		double num = 0;
	//		std::cout << "enter the number\n";
	//		std::cin >> num;
	//		second.SetElement(i + 1, j + 1, num);
	//	}
	//}
	//Matrix multiplication = temp.GetMatrixMultiplicationNew(second);
	//return 0;

	std::string pathOfFileToRead("C:\\Users\\digvijays\\Desktop\\GeometryAssignments\\GeometryAssignmentsModified\\GeometryAssignments\\PointsFromFile.txt");
	std::string pathOfFileToWrite("C:\\Users\\digvijays\\Desktop\\GeometryAssignments\\GeometryAssignmentsModified\\GeometryAssignments\\PointsToFile.txt");

	FileHandling fileHandling(pathOfFileToRead, pathOfFileToWrite);

	std::vector<Point3D> pointsFromFile;
	fileHandling.ReadPointsFromFile(pointsFromFile);

	int choice = 0;
	std::cout << "Enter your choce\n";
	std::cin >> choice;

	std::vector<Point3D> pointOnCurve;

	switch (choice)
	{
	case 1:  //Assignment 1: Get points along polyline and dump in Rhino                    
	{
		PolyLine3D polyLine(pointsFromFile);
		polyLine.GetPointsAlongSegments(pointOnCurve);
		
		break;
	}
	case 2:  //Assignment 2: Get approximate projection point on segment and dump in Rhino 
	{
		Point3D pointToProject(pointsFromFile[0]);
		Segment3D segment(pointsFromFile[1], pointsFromFile[2]);
		segment.GetPointsOnSegment(pointOnCurve);
		pointOnCurve.push_back(segment.GetApproximateProjectionOfPointOnSegment(pointToProject));
		break;
	}
	case 3:  //Assignment 3: Get precise projection point on segment and dump in Rhino 
	{
		Point3D pointToProject(pointsFromFile[0]);
		Segment3D segment(pointsFromFile[1], pointsFromFile[2]);
		pointOnCurve.push_back(segment.GetPreciseProjectionOfPointOnSegmentUsingVector(pointToProject));
		break;
	}
	case 4:  //Assignment 4: Get precise projection point along polyline and dump in Rhino 
	{
		PolyLine3D polyline(pointsFromFile);
		polyline.GetPointsAlongSegments(pointOnCurve);
		Point3D point(4,2,3);
		Point3D precisePoint =polyline.GetPrecisePointOfProjection(point);
		pointOnCurve.push_back(precisePoint);
		break;
	}
	case 5:  //Assignment 5: Get approx projection point along polyline and dump in Rhino 
	{
		PolyLine3D polyline(pointsFromFile);
		polyline.GetPointsAlongSegments(pointOnCurve);
		Point3D point(4, 2, 3);
		Point3D precisePoint = polyline.GetApproximateOfProjection(point);
		pointOnCurve.push_back(precisePoint);
		break;
	}
	case 6:   //Assignment 6: form a polycurve with the points given by user and dump points along the curve in Rhino 
	{
		Polynomial polynomial(pointsFromFile);
		polynomial.GetPointsAlongPolycurve(pointOnCurve);
		break;
	}
	case 7:   //Assignment 7: form a cubic Hermite Curve with the points and tangents given by user and dump points along the curve in Rhino 
	{
		BezierCurve bezierCurve(pointsFromFile);
		bezierCurve.GetPoinstAlongBazierCurve(1000, pointOnCurve);
		break;
	}
	case 8:    //Assignment 8: form a Bezier Curve with the control points given by user and dump points along the curve in Rhino 
	{
		std::vector<std::pair<Point3D, Vector3D>> temp;
		temp.push_back(std::make_pair(Point3D(1,0,0),Vector3D(0,1,0)));
		temp.push_back(std::make_pair(Point3D(2,1,1),Vector3D(1, 0, 0)));
		temp.push_back(std::make_pair(Point3D(3, 0, 2),Vector3D(0, -1, 0)));
		temp.push_back(std::make_pair(Point3D(4, 1, 3),Vector3D(-1, 0, 1)));
		CubicHermiteCurve hermiteCurve(temp);
		hermiteCurve.GetPointsAlongCubitHermiteCurve(pointOnCurve);
		break;
	}
	case 9:   //Assignment 7: form a BSpline Curve with given control points and degree given by user and dump points along the curve in Rhino 
	{
		int degree = 0;
		std::cout << "Enter degree for bspline curve\n";
		std::cin >> degree;
		BSplineCurve bsplineCurve(pointsFromFile, degree);
		bsplineCurve.GetPointsALongBSplineCurve(pointOnCurve);
		break;
	}
	default:
	{
		std::cout << "Invalid Choice!!!\n";
	}
	}

	fileHandling.WritePointsToFile(pointOnCurve);

	return 0;
}
