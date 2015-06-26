#ifndef POINT_UTIL_H
#define POINT_UTIL_H

#include <vector>
#include <set>
#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DepthFirstVisitor.h"

namespace PointUtil {
	template <typename Point>
	bool areAlmostSimilar(const Point& point, const Point& other);


	template <typename Domain, typename Point>
	Domain computeBoundingBox(const std::vector<Point>& points);


	/*
	 * We have to visit the direct neighbours in order to have a container with voxels
	 * ordered sequentially by their connexity
	 * Otherwise we have a point container with points which are not neighbours
	 * and this impairs maximal segment recognition
	 */
	template <typename Container, typename Image, typename Point>
	Container containerFromDepthTraversal(const Image& image, const Point& point, int thresholdMin,
		int thresholdMax);

	template <typename Point>
	std::vector<Point> linkTwoPoints(const Point& first, const Point& second);
}

template <typename Point>
bool PointUtil::areAlmostSimilar(const Point& point, const Point& other) {
	typename Point::Scalar otherx = other[0];
	typename Point::Scalar othery = other[1];
	typename Point::Scalar otherz = other[2];

	typename Point::Scalar pointx = point[0];
	typename Point::Scalar pointy = point[1];
	typename Point::Scalar pointz = point[2];
	
	bool sameX = pointx == otherx || pointx == otherx + 1 || pointx == otherx-1;
	bool sameY = pointy == othery || pointy == othery + 1 || pointy == othery-1;
	bool sameZ = pointz == otherz || pointz == otherz + 1 || pointz == otherz-1;
	return sameX && sameY && sameZ;
}

template < typename Domain, typename Point>
Domain PointUtil::computeBoundingBox(const std::vector<Point> & points) {
	int maximum = std::numeric_limits<int>::max();
	int min_x = maximum, min_y = maximum, min_z = maximum;
	int max_x = -maximum, max_y = -maximum, max_z = -maximum;
	for (const Point & point : points) {
		min_x = point[0] < min_x ? point[0] : min_x;
		min_y = point[1] < min_y ? point[1] : min_y;
		min_z = point[2] < min_z ? point[2] : min_z;
		max_x = point[0] > max_x ? point[0] : max_x;
		max_y = point[1] > max_y ? point[1] : max_y;
		max_z = point[2] > max_z ? point[2] : max_z;
	}
	Domain domain({min_x, min_y, min_z}, {max_x, max_y, max_z});
	return domain;
}

template <typename Container, typename Image, typename Point>
Container PointUtil::containerFromDepthTraversal(const Image& image, const Point& point, int thresholdMin,
										  int thresholdMax) {
	using namespace DGtal;

	
	typedef MetricAdjacency<Z3i::Space, 3> Graph;
	typedef DepthFirstVisitor<Graph, std::set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;
  
	Graph graph;
	Visitor visitor( graph, point );
	MyNode node;
	Container container;
    
	while ( !visitor.finished() ) 
	{
		node = visitor.current();
		if ( image.domain().isInside(node.first) &&
			 image(node.first) >= thresholdMin &&
			 image(node.first) <= thresholdMax ) { //is inside domain
		    container.push_back(node.first);
			visitor.expand();
		}
		else
			visitor.ignore();
	}
	return container;
}


template <typename Point>
std::vector<Point> PointUtil::linkTwoPoints(const Point& first, const Point& second) {

	Point source = first, destination = second;
	
	std::vector<Point> pointsBetween;
	DGtal::Z3i::RealPoint slope = second - first;
	
	std::set<int> coordinates = {0,1,2};
	double maxValue = 0;
	int indexForMaxValue = 0;
	for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
		double value = abs(slope[*it]);
		if (value > maxValue) {
			maxValue = value;
			indexForMaxValue = *it;
		}
	}

	std::set<int> otherCoordinates;
	for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
		if (*it != indexForMaxValue)
			otherCoordinates.insert(*it);
	}

	int otherIndex1 = *otherCoordinates.begin();
	int otherIndex2 = *(++otherCoordinates.begin());
	int valOtherValue1 = slope[otherIndex1] == 0 ? 0 : trunc(slope[indexForMaxValue]/slope[otherIndex1]);
	int valOtherValue2 = slope[otherIndex2] == 0 ? 0 : trunc(slope[indexForMaxValue]/slope[otherIndex2]);

	if (first[indexForMaxValue] > second[indexForMaxValue]) {
		source = second;
	 	destination = first;
		slope = -slope;
	}

	int currentValue1 = source[otherIndex1];
	int currentValue2 = source[otherIndex2];	
	
	for (int i = source[indexForMaxValue]+1; i < destination[indexForMaxValue]; i++) {
		Point p;
		p[indexForMaxValue] = i;
		p[otherIndex1] = currentValue1;
		p[otherIndex2] = currentValue2;
		pointsBetween.push_back(p);
		if (valOtherValue1 != 0 && i != 0 && i % valOtherValue1 == 0) {
			currentValue1 = (slope[otherIndex1] > 0 ? currentValue1+1 : currentValue1 - 1);
			Point pprime;
			pprime[indexForMaxValue] = i;
			pprime[otherIndex1] = currentValue1;
			pprime[otherIndex2] = currentValue2;
			pointsBetween.push_back(pprime);
		}
		if (valOtherValue2 != 0 && i != 0  && i % valOtherValue2 == 0) {
			currentValue2 = (slope[otherIndex2] > 0 ? currentValue2+1 : currentValue2 - 1);
			Point pprime;
			pprime[indexForMaxValue] = i;
			pprime[otherIndex1] = currentValue1;
			pprime[otherIndex2] = currentValue2;
			pointsBetween.push_back(pprime);
		}
	}
	return pointsBetween;
}

#endif
