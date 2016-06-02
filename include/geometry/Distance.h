#ifndef DISTANCE_H
#define DISTANCE_H

#include <math.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
 
inline double euclideanDistance(float x1, float y1, float x2, float y2) {
	return sqrt( pow( (x2 - x1), 2) + pow( (y2 - y1), 2));
}

template <typename Point>
inline double euclideanDistance(Point p, Point other) {
	return sqrt(pow( (p[0] - other[0]), 2) + pow( (p[1] - other[1]), 2) + pow( (p[2] - other[2]), 2) );
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
#endif
