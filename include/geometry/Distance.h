#ifndef DISTANCE_H
#define DISTANCE_H

#include <math.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

namespace Distance {
	inline double euclideanDistance(float x1, float y1, float x2, float y2) {
		return sqrt( pow( (x2 - x1), 2) + pow( (y2 - y1), 2));
	}

	template <typename Point>
	inline double euclideanDistance(Point p, Point other) {
		int dim = Point::dimension;
		double sum = 0;
		for (int i = 0; i < dim; i++) {
			sum += pow( (p[i] - other[i]), 2);
		}
		return sqrt(sum);
	}

	template <typename Container>
	inline double hausdorffDistance(const Container& first, const Container& second) {
		double distanceMax = 0;
		for (auto it = first.begin(), ite = first.end(); it != ite; ++it) {
			DGtal::Z3i::Point firstP = *it;
			DGtal::Z3i::Point closestPointInTheoretical =  *min_element(second.begin(), second.end(), [&](const DGtal::Z3i::Point& one, const DGtal::Z3i::Point& two) {
					return DGtal::Z3i::l2Metric(one, firstP) < DGtal::Z3i::l2Metric(two, firstP);
				});
			double distance = DGtal::Z3i::l2Metric(closestPointInTheoretical, firstP);
			if (distance > distanceMax)
				distanceMax = distance;
		}
		return distanceMax;
	}

	template <typename Container>
	inline double distanceSet(const Container& first, const Container& second) {
		double distanceMax = std::numeric_limits<double>::max();
		for (auto it = first.begin(), ite = first.end(); it != ite; ++it) {
			DGtal::Z3i::Point firstP = *it;
			DGtal::Z3i::Point closestPointInTheoretical =  *min_element(second.begin(), second.end(), [&](const DGtal::Z3i::Point& one, const DGtal::Z3i::Point& two) {
					return DGtal::Z3i::l2Metric(one, firstP) < DGtal::Z3i::l2Metric(two, firstP);
				});
			double distance = DGtal::Z3i::l2Metric(closestPointInTheoretical, firstP);
			if (distance < distanceMax)
				distanceMax = distance;
		}
		return distanceMax;
	}

	template <typename Point, typename Container>
	inline double geodesicDistance(const Point& first, const Point& second, const Container& object) {
		typedef DGtal::BreadthFirstVisitor<DGtal::Z3i::Object26_6, std::set<Point> > Visitor;
		typedef typename Visitor::Node MyNode;
		DGtal::Z3i::Object26_6 graph(DGtal::Z3i::dt26_6, object);
		Visitor visitor( graph, first );
		MyNode node;
		double distance = std::numeric_limits<double>::max();
		while ( !visitor.finished() )
		{
			node = visitor.current();
		    if (node.first == second) {
				distance = node.second;
				break;
			}
			visitor.expand();
		}
		return distance;
	}
};
#endif
