#include "FileHandling.h"
#include "Point3D.h"
#include "Segment3D.h"
#include "Vector3D.h"
#include <fstream>
#include <sstream>

//-----------------------------------------------------------------------------

FileHandling::FileHandling(const std::string& filePathToReadPoints, 
	                       const std::string& filePathToWritePoits)
{
	m_filePathToReadPoints  = filePathToReadPoints;
	m_filePathToWritePoints = filePathToWritePoits;
}

//-----------------------------------------------------------------------------

void FileHandling::ReadPointsFromFile(std::vector<Point3D>& pointsFromFile)const
{
	std::fstream readFile(m_filePathToReadPoints,std::ios::in);
	if (readFile.is_open())
	{
		std::string currentLine;
		while (std::getline(readFile,currentLine))
		{
			std::istringstream tempLine(currentLine);
			std::vector<std::string> point;
			std::string currentCoordinate;
			while (std::getline(tempLine,currentCoordinate,','))
			{
				point.push_back(currentCoordinate);
			}
			Point3D temp(stod(point[0]), stod(point[1]), stod(point[2]));
			pointsFromFile.push_back(temp);
		}
	}
	readFile.close();
}

//-----------------------------------------------------------------------------

void FileHandling::WritePointsToFile(const std::vector<Point3D>& pointsOnSegments)const
{
	std::fstream writeFile(m_filePathToWritePoints,std::ios::trunc|std::ios::out);
	if (writeFile.is_open())
	{
		for (int iPoint = 0; iPoint < pointsOnSegments.size(); ++iPoint)
		{
			Point3D currentPoint = pointsOnSegments[iPoint];
			writeFile << currentPoint.GetX() << ',' << currentPoint.GetY() << ',' << currentPoint.GetZ() << '\n';
		}
	}
	writeFile.close();
}

//-----------------------------------------------------------------------------

FileHandling::~FileHandling()
{

}

//-----------------------------------------------------------------------------
