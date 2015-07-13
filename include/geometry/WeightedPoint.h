#ifndef WEIGHTED_POINT_H
#define WEIGHTED_POINT_H

template <typename Point>
class WeightedPoint {
public:
	WeightedPoint() : myPoint(), myWeight(0.) {}
	WeightedPoint( const Point& aPoint, double aWeight ) : myPoint(aPoint), myWeight(aWeight) {}
	friend bool operator<(const WeightedPoint& it, const WeightedPoint& other) {
		if (it.myPoint == other.myPoint)
		    return false;
		if (it.myWeight >= other.myWeight)
			return true;
		return false;
	}

	friend bool operator!=(const WeightedPoint& it, const WeightedPoint& other) {
		return (it.myPoint != other.myPoint || it.myWeight != other.myWeight);
	}
public:
	Point myPoint;
	double myWeight;
};

#endif
