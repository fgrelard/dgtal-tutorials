#ifndef WEIGHTED_POINT_COUNT_H
#define WEIGHTED_POINT_COUNT_H

#include "geometry/WeightedPoint.h"

template <typename Point>
class WeightedPointCount : public WeightedPoint<Point> {
public:
	WeightedPointCount(const Point& aPoint) : WeightedPoint<Point>(aPoint, 0) {}
	double getWeightedValue() { return WeightedPoint<Point>::myWeight / myCount; }
public:
	int myCount = 0;
};

#endif
