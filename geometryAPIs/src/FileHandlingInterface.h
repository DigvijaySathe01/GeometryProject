#pragma once

#include <vector>
#include <string>

class Segment3D;
class Point3D;

class FileHandlingInterface
{

public:

	//Function to read points from file
	virtual void ReadPointsFromFile(std::vector<Point3D>& pointsFromFile)const = 0;

	//Function to write points on polyline to the file
	virtual void WritePointsToFile(const std::vector<Point3D>& pointsOnSegments)const = 0;

	//Destructor
	virtual ~FileHandlingInterface() =default;

};
