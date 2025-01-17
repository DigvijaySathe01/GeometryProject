#include "PolyLine3D.h"
#include "Point3D.h"
#include "Segment3D.h"
#include <limits>

//-----------------------------------------------------------------------------

PolyLine3D::PolyLine3D(const std::vector<Point3D>& pointsFromFile)
{
	Point3D prevPoint = pointsFromFile[0];
	for (int iPoint = 1; iPoint < pointsFromFile.size(); ++iPoint)
	{
		Point3D currentPoint = pointsFromFile[iPoint];
		Segment3D* segment = new Segment3D(prevPoint,currentPoint);
		m_segments.push_back(segment);
		prevPoint = currentPoint;
	}
}

//-----------------------------------------------------------------------------

void PolyLine3D::GetPointsAlongSegment(std::vector<Point3D>& pointsOnSegments, int index)
{
	if (index < pointsOnSegments.size())
		m_segments[index]->GetPointsOnSegment(pointsOnSegments);
}

//-----------------------------------------------------------------------------

void PolyLine3D::GetPointsAlongSegments(std::vector<Point3D>& pointsOnSegments)
{
	for (int iPoint = 0; iPoint < m_segments.size(); ++iPoint)
	{
		m_segments[iPoint]->GetPointsOnSegment(pointsOnSegments);
	}
}

//-----------------------------------------------------------------------------

Point3D PolyLine3D::GetPrecisePointOfProjection(const Point3D& pointToProject) const
{
	std::vector<Point3D> precisePointsOnSegments;
	for (int iSegment = 0; iSegment < m_segments.size(); ++iSegment)
	{
		Point3D pointOfProjection;
		pointOfProjection = m_segments[iSegment]->GetPreciseProjectionOfPointOnSegmentUsingVector(pointOfProjection);
		precisePointsOnSegments.push_back(pointOfProjection);
	}

	double distance = std::numeric_limits<double>::max();
	Point3D precisePointOfProjection;
	for (int iPoint = 0; iPoint < precisePointsOnSegments.size(); ++iPoint)
	{
		if (pointToProject.GetDistance(precisePointsOnSegments[iPoint]) <= distance)
		{
			precisePointOfProjection = precisePointsOnSegments[iPoint];
			distance = pointToProject.GetDistance(precisePointsOnSegments[iPoint]);
		}
	}
	return precisePointOfProjection;
}

//-----------------------------------------------------------------------------

Point3D PolyLine3D::GetApproximateOfProjection(const Point3D& pointToProject) const
{
	std::vector<Point3D> precisePointsOnSegments;
	for (int iSegment = 0; iSegment < m_segments.size(); ++iSegment)
	{
		Point3D pointOfProjection;
		pointOfProjection = m_segments[iSegment]->GetApproximateProjectionOfPointOnSegment(pointOfProjection);
		precisePointsOnSegments.push_back(pointOfProjection);
	}

	double distance = std::numeric_limits<double>::max();
	Point3D precisePointOfProjection;
	for (int iPoint = 0; iPoint < precisePointsOnSegments.size(); ++iPoint)
	{
		if (pointToProject.GetDistance(precisePointsOnSegments[iPoint]) <= distance)
		{
			precisePointOfProjection = precisePointsOnSegments[iPoint];
			distance = pointToProject.GetDistance(precisePointsOnSegments[iPoint]);
		}
	}
	return precisePointOfProjection;
}

//-----------------------------------------------------------------------------


PolyLine3D::~PolyLine3D()
{
	for (int iSize = 0; iSize < m_segments.size(); ++iSize)
	{
		delete m_segments[iSize];
	}
}

//-----------------------------------------------------------------------------