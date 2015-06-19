#ifndef WEIGHTED_POINT_H
#define WEIGHTED_POINT_H

template <typename Point>
class WeightedPoint {
public:
	WeightedPoint( const Point& aPoint, double aWeight ) : myPoint(aPoint), myWeight(aWeight) {}
	friend bool operator<(const WeightedPoint& it, const WeightedPoint& other) {
		if (it.d >= other.d) {
			return true;
		}
		return false;
	}
public:
	Point myPoint;
	double myWeight;
};

#endif
