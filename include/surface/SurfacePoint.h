#ifndef __SURFACE_POINT__
#define __SURFACE_POINT__

#include <vector>
#include <algorithm>

template <typename Point, typename Surfel>
class SurfacePoint : public Point {
public:
	typedef typename Point::Component Scalar;
public:
	SurfacePoint() : Point(), myNumberOfSurfels{0} {}
	SurfacePoint(const Scalar& x, const Scalar& y, const Scalar& z) : Point{x,y,z},myNumberOfSurfels{0} {}
	SurfacePoint(const Point& point, int numberOfSurfels, const std::vector<Surfel>& surfels) : Point{point[0], point[1], point[2]}, myNumberOfSurfels{numberOfSurfels}, mySurfels{surfels} {}
	bool isScanned() { return myVisitedSurfels >= myNumberOfSurfels;}
	void oneSurfelVisit() { myVisitedSurfels++; }
	Surfel& surfelToScan() { return mySurfels[myVisitedSurfels]; }
	int getNumberOfSurfels() { return myNumberOfSurfels; }
	bool contains(const Surfel& s) {
		return (std::find(mySurfels.begin(), mySurfels.end(), s) != mySurfels.end());
	}

	std::vector<Surfel> surfelsAtPoint() { return mySurfels; }
private:
    int myVisitedSurfels = 1;
	int myNumberOfSurfels;
	std::vector<Surfel> mySurfels;
};

#endif
