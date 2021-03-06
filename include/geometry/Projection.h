#ifndef PROJECTION_H
#define PROJECTION_H

#include <map>
#include <vector>
#include "DGtal/helpers/StdDefs.h"

namespace Projection {
	/**
	 * @param referencePoint the point to test
	 * @param points the set of points to project onto
	 * @return for each point in points the associated euclidean distance to referencePoint
	 */
	template <typename Point>
	std::map<Point, double> constructDistanceMapFromReferencePoint(const Point& referencePoint, const DGtal::Z3i::DigitalSet& points);
	
	/** 
	 * @param referencePoint the point to test
	 * @param points the set of points to project onto
	 * @return a list of points that correspond to the projections of referencePoint on points
	 */
	template <typename Point>
	DGtal::Z3i::DigitalSet extractFeaturePoints(const Point& referencePoint, const DGtal::Z3i::DigitalSet& points);
}

template <typename Point>
std::map<Point, double> Projection::constructDistanceMapFromReferencePoint(const Point& referencePoint, const DGtal::Z3i::DigitalSet& points) {
	std::map<Point, double> distanceMap;

	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		distanceMap[*it] = DGtal::Z3i::l2Metric(referencePoint, *it);
	}

	return distanceMap;
}

template <typename Point>
DGtal::Z3i::DigitalSet Projection::extractFeaturePoints(const Point& referencePoint, const DGtal::Z3i::DigitalSet& points)
{
	if (points.find(referencePoint) != points.end()) return DGtal::Z3i::DigitalSet(points.domain());
	//Computing neighborhood around referencePoint
	typedef DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> MetricAdjacency;
	std::vector<Point> neighbours26;
	std::back_insert_iterator<std::vector<Point>> outIt(neighbours26);
	MetricAdjacency::writeNeighbors(outIt, referencePoint);

	//Extracting feature points
	DGtal::Z3i::DigitalSet featurePoints(points.domain());

	std::map<Point, double> distanceMap = Projection::constructDistanceMapFromReferencePoint(referencePoint, points);
	
	double cmap_min = min_element(distanceMap.begin(), distanceMap.end(), [&](const std::pair<Point, double>& one, const std::pair<Point, double>& two) {
			return one.second < two.second;
		})->second;
	
	for (auto it = distanceMap.begin(), ite = distanceMap.end(); it != ite; ++it) {
		if (it->second == cmap_min)
			featurePoints.insert(it->first);
	}
	for (auto it = neighbours26.begin(), ite = neighbours26.end(); it != ite; ++it) {
		std::map<Point, double> distanceMap = Projection::constructDistanceMapFromReferencePoint(*it, points);
		
		double cmap_min = min_element(distanceMap.begin(), distanceMap.end(), [&](const std::pair<Point, double>& one, const std::pair<Point, double>& two) {
			return one.second < two.second;
			})->second;
		
		for (auto it = distanceMap.begin(), ite = distanceMap.end(); it != ite; ++it) {
			if (it->second == cmap_min)
				featurePoints.insert(it->first);
		}
	}
	return featurePoints;
}



#endif
