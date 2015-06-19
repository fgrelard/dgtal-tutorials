#ifndef DISTANCE_H
#define DISTANCE_H

#include <math.h>

inline double euclideanDistance(float x1, float y1, float x2, float y2) {
	return sqrt( pow( (x2 - x1), 2) + pow( (y2 - y1), 2));
}

template <typename Point>
inline double euclideanDistance(Point p, Point other) {
	return sqrt(pow( (p[0] - other[0]), 2) + pow( (p[1] - other[1]), 2) + pow( (p[2] - other[2]), 2) );
}

#endif
