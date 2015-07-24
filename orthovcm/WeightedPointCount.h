#ifndef WEIGHTED_POINT_COUNT_H
#define WEIGHTED_POINT_COUNT_H

#include "geometry/WeightedPoint.h"

template <typename Point>
class WeightedPointCount : public WeightedPoint<Point> {
	typedef WeightedPoint<Point> Base;
public:
	using Base::Base;
	double getWeightedValue() { return WeightedPoint<Point>::myWeight / myCount; }
	friend bool operator<(const WeightedPointCount& it, const WeightedPointCount& other) {
		return (it.myWeight >= other.myWeight);
	}
public:
	int myCount = 0;
};

#endif
