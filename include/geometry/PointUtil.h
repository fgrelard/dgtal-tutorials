#ifndef POINT_UTIL_H
#define POINT_UTIL_H

#include <vector>
#include <limits>

namespace PointUtil {
	template <typename Point>
	bool areAlmostSimilar(const Point& point, const Point& other);


	template <typename Domain, typename Point>
	Domain computeBoundingBox(const std::vector<Point>& points);
}

template <typename Point>
bool PointUtil::areAlmostSimilar(const Point& point, const Point& other) {
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

template < typename Domain, typename Point>
Domain PointUtil::computeBoundingBox(const std::vector<Point> & points) {
	int maximum = std::numeric_limits<int>::max();
	int min_x = maximum, min_y = maximum, min_z = maximum;
	int max_x = -maximum, max_y = -maximum, max_z = -maximum;
	for (const Point & point : points) {
		min_x = point[0] < min_x ? point[0] : min_x;
		min_y = point[1] < min_y ? point[1] : min_y;
		min_z = point[2] < min_z ? point[2] : min_z;
		max_x = point[0] > max_x ? point[0] : max_x;
		max_y = point[1] > max_y ? point[1] : max_y;
		max_z = point[2] > max_z ? point[2] : max_z;
	}
	Domain domain({min_x, min_y, min_z}, {max_x, max_y, max_z});
	return domain;
}


#endif
