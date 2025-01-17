#pragma once
#include <vector>
#include <string>

#include "FileHandlingInterface.h"
class Segment3D;
class Point3D;

class FileHandling : public FileHandlingInterface
{
private:
	std::string m_filePathToReadPoints;       //path of the file which contains points to form polyline
	std::string m_filePathToWritePoints;      //path of the file to which we will write the points

public:

	//Parameterized constructor
	FileHandling(const std::string& filePathToReadPoints, const std::string& filePathToWritePoits);

	//Function to read points from file
	void ReadPointsFromFile(std::vector<Point3D>& pointsFromFile)const override;

	//Function to write points on polyline to the file
	void WritePointsToFile(const std::vector<Point3D>& pointsOnSegments)const override;

	//Destructor
	~FileHandling();
};

