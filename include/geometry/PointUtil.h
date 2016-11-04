#ifndef POINT_UTIL_H
#define POINT_UTIL_H

#include <vector>
#include <set>
#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "Polygon.h"

namespace PointUtil {
	template <typename Point>
	bool areNeighbors(const Point& point, const Point& other);


	template <typename Domain, typename Container>
	Domain computeBoundingBox(const Container& points);


	/*
	 * We have to visit the direct neighbours in order to have a container with voxels
	 * ordered sequentially by their connexity
	 * Otherwise we have a point container with points which are not neighbours
	 * and this impairs maximal segment recognition
	 */
	template <typename Container, typename Image, typename Point>
	Container containerFromDepthTraversal(const Image& image, const Point& point, int thresholdMin,
		int thresholdMax);

	template <typename Container, typename OtherContainer, typename Point>
	Container containerFromDepthTraversal(const OtherContainer& image, const Point& point);

	template <typename Point>
	std::vector<Point> linkTwoPoints(const Point& first, const Point& second);

	template <typename Point>
	std::vector<Point> linkTwoPoints26(const Point& first, const Point& second);

	template <typename Point>
	std::vector<Point> bezierCurve(const Point& source,
								   const Point& destination,
								   const Point& controlPoint1,
								   const Point& controlPoint2);


	template <typename Point>
	std::vector<Polygon<Point> > createPolygonSubdivision(const Polygon<Point>& polygon);

	template <typename Point>
	void recursivePolygonSubdivisionBezier(const Polygon<Point>& polygon, std::map<Point, Point>& points);

	template <typename Point>
	bool stopRecursivePolygonDecomposition(const Polygon<Point>& polygon);

	int computeDiameter(const std::set<int>& equation);

	template <typename Point>
	std::vector<std::pair<Point, Point> > orderPoints(const Point& source,
								   const Point& destination,
								   const std::map<Point, Point>& points);

	template <typename Point>
	std::vector<Point> bezierCurveDeCasteljau(const Point& source,
											  const Point& destination,
											  const Point& controlPoint1,
											  const Point& controlPoint2);
	
	template <typename Point, typename Container, typename Vector>
	Point trackPoint(const Point& initial, const Container& container, const Vector& vector);
}

template <typename Point>
bool PointUtil::areNeighbors(const Point& point, const Point& other) {
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

template <typename Domain, typename Container>
Domain PointUtil::computeBoundingBox(const Container & points) {
	int maximum = std::numeric_limits<int>::max();
	int min_x = maximum, min_y = maximum, min_z = maximum;
	int max_x = -maximum, max_y = -maximum, max_z = -maximum;
	for (const auto & point : points) {
		min_x = point[0] < min_x ? point[0] : min_x;
		min_y = point[1] < min_y ? point[1] : min_y;
		min_z = point[2] < min_z ? point[2] : min_z;
		max_x = point[0] > max_x ? point[0] : max_x;
		max_y = point[1] > max_y ? point[1] : max_y;
		max_z = point[2] > max_z ? point[2] : max_z;
	}
	Domain domain({min_x-1, min_y-1, min_z-1}, {max_x+1, max_y+1, max_z+1});
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

template <typename Container, typename OtherContainer, typename Point>
Container PointUtil::containerFromDepthTraversal(const OtherContainer& image, const Point& point) {
	using namespace DGtal;

	
	typedef DGtal::Z3i::Object26_6 Graph;
	typedef DepthFirstVisitor<Graph, std::set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;
  
	Graph graph(DGtal::Z3i::dt26_6, image);
	Visitor visitor( graph, point );
	MyNode node;
	Container container;
    
	while ( !visitor.finished() ) 
	{
		node = visitor.current();
		container.push_back(node.first);
		visitor.expand();
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
	float valOtherValue1 = slope[otherIndex1] == 0 ? 0 : fabs(1.0*slope[indexForMaxValue]/slope[otherIndex1]);
	float valOtherValue2 = slope[otherIndex2] == 0 ? 0 : fabs(1.0*slope[indexForMaxValue]/slope[otherIndex2]);

	bool toReverse = false;
	if (first[indexForMaxValue] > second[indexForMaxValue]) {
		source = second;
	 	destination = first;
		slope = -slope;
		toReverse = true;
	}

	int currentValue1 = source[otherIndex1];
	int currentValue2 = source[otherIndex2];	
	float errorValue1 = 1.0;
	float errorValue2 = 1.0;
	for (int i = source[indexForMaxValue]+1; i < destination[indexForMaxValue]; i++) {
		Point p;
		p[indexForMaxValue] = i;
		p[otherIndex1] = currentValue1;
		p[otherIndex2] = currentValue2;
		pointsBetween.push_back(p);
		if (valOtherValue1 != 0 && i != 0 && errorValue1 >= valOtherValue1) {
			currentValue1 = (slope[otherIndex1] > 0 ? currentValue1+1 : currentValue1 - 1);
			Point pprime;
			pprime[indexForMaxValue] = i;
			pprime[otherIndex1] = currentValue1;
			pprime[otherIndex2] = currentValue2;
			pointsBetween.push_back(pprime);
			errorValue1 -= (valOtherValue1 - 1);
		} else {
			errorValue1 += 1.0;
		}
		if (valOtherValue2 != 0 && i != 0  && errorValue2 >= valOtherValue2) {
			currentValue2 = (slope[otherIndex2] > 0 ? currentValue2+1 : currentValue2 - 1);
			Point pprime;
			pprime[indexForMaxValue] = i;
			pprime[otherIndex1] = currentValue1;
			pprime[otherIndex2] = currentValue2;
			pointsBetween.push_back(pprime);
			errorValue2 -= (valOtherValue2 - 1);
		} else
			errorValue2 += 1.0;
	}
	if (toReverse)
		std::reverse(pointsBetween.begin(), pointsBetween.end());
	return pointsBetween;
}

template <typename Point>
std::vector<Point> PointUtil::bezierCurve(const Point& source,
							   const Point& destination,
							   const Point& controlPoint1,
							   const Point& controlPoint2) {
	std::vector<Point> bezier;
	for (double t =0; t < 1; t+=0.01) {
		Point current;
		for (int k = 0; k < Point::dimension; k++) {
			double coordinateAtK = pow((1-t), 3) * source[k] +
				3 * pow((1-t), 2) * t * controlPoint1[k] +
				3 * (1-t) * pow(t, 2) * controlPoint2[k] +
				pow(t, 3) * destination[k];
			current[k] = std::round(coordinateAtK);
		}
		bezier.push_back(current);
	}
	return bezier;
}

template <typename Point>
std::vector<Polygon<Point> > PointUtil::createPolygonSubdivision(const Polygon<Point>& polygon) {
	typedef std::pair<Point, Point> Segment;

	Point source = polygon.getPolygon()[0];
	Point controlPoint1 = polygon.getPolygon()[1];
	Point controlPoint2 = polygon.getPolygon()[2];
	Point destination = polygon.getPolygon()[3];

	std::vector<Polygon<Point> > polygons;
    Segment s1(source, controlPoint1);
	Segment s2(controlPoint1, controlPoint2);
	Segment s3(controlPoint2, destination);

	std::vector<Segment> segments = {s1, s2, s3};
	std::vector<Point> newPoints;
	while (segments.size() > 0) {
		std::vector<Point> currentPoints;
		for (const Segment& s : segments) {
			Point first = s.first;
			Point second = s.second;
		    Point middle = (first + second) * 1/2.;
			newPoints.push_back(middle);
			currentPoints.push_back(middle);
		}

		segments.clear();
		for (int i = 0; i < currentPoints.size() - 1; i++) {
			Segment s(currentPoints[i], currentPoints[i+1]);
			segments.push_back(s);
		}
	}
	Polygon<Point> polygon1{source, newPoints[0], newPoints[3], newPoints[5]};
	Polygon<Point> polygon2{newPoints[5], newPoints[4], newPoints[2], destination};
	polygons.push_back(polygon1);
	polygons.push_back(polygon2);
	return polygons;
}



int PointUtil::computeDiameter(const std::set<int>& equation) {
	return *equation.rbegin() - *equation.begin();
}

template <typename Point>
bool PointUtil::stopRecursivePolygonDecomposition(const Polygon<Point>& polygon) {
	std::vector<Point> points = polygon.getPolygon();
	Point dirVector = points[3] -points[0];
	std::set<int> firstEquation, secondEquation;
	for (const Point& point : points) {
		int value = -dirVector[2] * point[0] + dirVector[0] * point[2];
		int secondValue = -dirVector[2] * point[1] + dirVector[1] * point[2];
	    firstEquation.insert(value);
		secondEquation.insert(secondValue);
	}
	double diameter1 = computeDiameter(firstEquation);
	double diameter2 = computeDiameter(secondEquation);

	double valueToCheck = std::max(std::abs(dirVector[0]), std::max(std::abs(dirVector[1]), std::abs(dirVector[2])));
	valueToCheck *= 4 / 3.;
	return (diameter1 <= valueToCheck && diameter2 <= valueToCheck);
}

template <typename Point>
void PointUtil::recursivePolygonSubdivisionBezier(const Polygon<Point>& polygon, std::map<Point, Point>& points) {
	if (!stopRecursivePolygonDecomposition(polygon)) {
	    std::vector<Polygon<Point> > polygons = createPolygonSubdivision(polygon);
		recursivePolygonSubdivisionBezier<Point>(polygons[0], points);
		recursivePolygonSubdivisionBezier<Point>(polygons[1], points);
	} else {
		points[polygon.getPolygon()[0]] = polygon.getPolygon()[3];
	}
}

template <typename Point>
std::vector<std::pair<Point, Point> > PointUtil::orderPoints(const Point& source,
											  const Point& destination,
											  const std::map<Point, Point>& points) {
	Point next = source;
	std::vector<std::pair<Point, Point>> orderedPoints;
	while (next != destination) {
		Point current = next;
		auto iterator = points.find(next);
		if (iterator == points.end()) break;
		next = iterator->second;
		orderedPoints.push_back(make_pair(current,  next));
	}
	return orderedPoints;
}


template <typename Point>
std::vector<Point> PointUtil::bezierCurveDeCasteljau(const Point& source,
													 const Point& destination,
													 const Point& controlPoint1,
													 const Point& controlPoint2) {
	Polygon<Point> polygon{source, controlPoint1, controlPoint2, destination};
	std::map<Point, Point> points;
	recursivePolygonSubdivisionBezier(polygon, points);
	std::vector<std::pair<Point, Point> > orderedPoints = PointUtil::orderPoints(source, destination, points);
	std::vector<Point> linkPoints;
	for (const auto & pair : orderedPoints) {
		std::vector<Point> link = linkTwoPoints(pair.first, pair.second);
		if (pair.first != source && pair.first != destination)
			linkPoints.push_back(pair.first);
		linkPoints.insert(linkPoints.end(), link.begin(), link.end());
		if (pair.second != destination && pair.second != source)
			linkPoints.push_back(pair.second);
	}
	return linkPoints;
}

template <typename Point>
std::vector<Point> PointUtil::linkTwoPoints26(const Point& first, const Point& second) {
	Point source = first, destination = second;

	std::vector<Point> pointsBetween;
	DGtal::Z3i::RealPoint slope = second - first;
	std::set<int> coordinates = {0,1,2};
	double maxValue = 0;
	int indexForMaxValue = 0, indexForSecond = 0, indexForThird = 0;
	for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
		double value = std::abs(slope[*it]);
		if (value > maxValue) {
			maxValue = value;
			indexForMaxValue = *it;
		}
	}

	maxValue = -1;
	for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
		if (*it == indexForMaxValue) continue;
		double value = std::abs(slope[*it]);
		if (value > maxValue) {
			maxValue = value;
			indexForSecond = *it;
		}
	}

	for (auto it = coordinates.begin(), ite = coordinates.end(); it != ite; ++it) {
		if (*it == indexForMaxValue || *it == indexForSecond) continue;
		indexForThird = *it;
	}


	bool toReverse = false;
	if (first[indexForMaxValue] > second[indexForMaxValue]) {
		source = second;
	 	destination = first;
		slope = -slope;
		toReverse = true;
	}

	//y-init
	double deltaX = slope[indexForMaxValue];
	int signX = (slope[indexForMaxValue] > 0) - (slope[indexForMaxValue] < 0);

	double deltaY = std::abs(slope[indexForSecond]);
	int signY = (slope[indexForSecond] > 0) - (slope[indexForSecond] < 0);
	double ey = 2 * deltaY - deltaX;
	double yinc1 = 2 * deltaY;
	double yinc2 = 2 * (deltaY - deltaX);

	//z-init
	double deltaZ = std::abs(slope[indexForThird]);
	int signZ = (slope[indexForThird] > 0) - (slope[indexForThird] < 0);
	double ez = 2 * deltaZ - deltaX;
	double zinc1 = 2 * deltaZ;
	double zinc2 = 2 * (deltaZ - deltaX);


	Point current = source, next = source;

	Point offsetX, offsetY, offsetZ;
	offsetX[indexForMaxValue] = signX;
	offsetY[indexForSecond] = signY;
	offsetZ[indexForThird] = signZ;

	int i = current[indexForMaxValue];
	int imax = destination[indexForMaxValue];
	while (i < imax) {
		i++;
		if (current != source && current != destination)
			pointsBetween.push_back(current);
		current += offsetX;
		if (ey < 0) {
			ey += yinc1;

		}
		else {
			ey += yinc2;
			current += offsetY;
		}

		if (ez < 0) {
			ez += zinc1;
		}
		else {
			ez += zinc2;
			current += offsetZ;
		}
	}
	if (toReverse)
		std::reverse(pointsBetween.begin(), pointsBetween.end());
	return pointsBetween;
}

template <typename Point, typename Container, typename Vector>
Point PointUtil::trackPoint(const Point& initial, const Container& container, const Vector& vector) {
	Point point(initial);
	int scalar = 1;
	while (container.find(point) != container.end()) {
		point = initial + vector*scalar;
		scalar++;
	}
	return point;
}

#endif
