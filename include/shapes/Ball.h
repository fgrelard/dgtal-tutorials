#ifndef BALL_H
#define BALL_H

#include <vector>
#include <set>
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "geometry/Distance.h"

template <typename Point>
class Ball {
public:
	Ball() : myRadius(0.0), myCenter({0,0,0}) {}
	Ball(const Point& center, double radius) : myCenter(center), myRadius(radius) {}
	Ball(const Ball& other) : myCenter(other.myCenter), myRadius(other.myRadius) {}
	bool contains(const Point& point) const {return euclideanDistance(point, myCenter) <= myRadius;}
	std::vector<Point> surfaceIntersection(const DGtal::Z3i::DigitalSet& setSurface);
	std::vector<Point> pointsInBall() const;
	std::vector<Point> pointsInHalfBall() const;
	DGtal::Z3i::DigitalSet pointsInBallSet() const;
	
	template <typename RealPoint>
	std::vector<Point> pointsInHalfBall(const RealPoint& normal) const;

	std::set<Point> pointsSurfaceBall() const;
	bool operator!=(const Ball & other) const {return (myCenter != other.myCenter || myRadius != other.myRadius);}
	Point getCenter()  {return myCenter;}
private:
	Point myCenter;
	double myRadius;
};

template <typename Point>
std::vector<Point> Ball<Point>::surfaceIntersection(const DGtal::Z3i::DigitalSet& setSurface) {
	std::vector<Point> intersection;
	for (auto it = setSurface.begin(), ite = setSurface.end(); it != ite; ++it) {
		double distance = euclideanDistance(*it, myCenter);
		if (distance >= myRadius-1 && distance <= myRadius) {
			intersection.push_back(*it);
		}
	}
	return intersection;
}

template <typename Point>
DGtal::Z3i::DigitalSet Ball<Point>::pointsInBallSet() const {
    Point lower(-myRadius + myCenter[0], -myRadius + myCenter[1], -myRadius + myCenter[2]);
	Point upper(myRadius + myCenter[0] + 1, myRadius + myCenter[1] + 1, myRadius + myCenter[2] + 1);
	DGtal::Z3i::DigitalSet::Domain domain(lower, upper);
	DGtal::Z3i::DigitalSet points(domain);
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);
				if (contains(p)) {
					points.insert(p);
				}
			}
		}
	}
	return points;
}

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
template <typename RealPoint>
std::vector<Point> Ball<Point>::pointsInHalfBall(const RealPoint& normal) const {
	std::vector<Point> points;
	double d = myCenter[0] * normal[0] + myCenter[1] * normal[1] + myCenter[2] * normal[2];
	for (typename Point::Scalar i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
		for (typename Point::Scalar j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
			for (typename Point::Scalar k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
				Point p(i, j, k);				
				double eq = p[0] * normal[0] + p[1] * normal[1] + p[2] * normal[2] - d;
				if (contains(p) && eq > 0) {
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
