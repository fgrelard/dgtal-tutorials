#ifndef BALL_H
#define BALL_H

#include <vector>
#include <set>
#include <iostream>

template <typename Point>
class Ball {
public:
	Ball() : myRadius(0.0), myCenter({0,0,0}) {}
	Ball(const Point& center, double radius) : myCenter(center), myRadius(radius) {}
	bool contains(const Point& point) const {return euclideanDistance(point, myCenter) <= myRadius;}
	std::vector<Point> pointsInBall() const;
	std::vector<Point> pointsInHalfBall() const;
	std::set<Point> pointsSurfaceBall() const;
	bool operator!=(const Ball & other) const {return (myCenter != other.myCenter || myRadius != other.myRadius);}
	Point getCenter()  {return myCenter;}
private:
	Point myCenter;
	double myRadius;
};

template <typename Point>
std::vector<Point> Ball<Point>::pointsInBall() const {	
	std::vector<Point> points;
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				if (contains(p)) {
					points.push_back(p);
				}
			}
		}
	}
	return points;	
}

template <typename Point>
std::vector<Point> Ball<Point>::pointsInHalfBall() const {	
	std::vector<Point> points;
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend =myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				if (contains(p)) {
					points.push_back(p);
				}
			}
		}
	}
	return points;	
}

template <typename Point>
std::set<Point> Ball<Point>::pointsSurfaceBall() const {
	std::set<Point> points;
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				double distance = euclideanDistance(p, myCenter);
				if (distance >= myRadius-1 && distance <= myRadius) {
					points.insert(p);
				}
			}
		}
	}
	return points;
}
#endif
