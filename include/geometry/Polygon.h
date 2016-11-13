#ifndef POLYGON_H
#define POLYGON_H

#include <initializer_list>
#include <vector>

template <typename Point>
class Polygon {
public:
        Polygon() { }
        Polygon(std::initializer_list<Point> l) {
                myPoints = l;
        }

        std::vector<Point> getPolygon() const { return myPoints; }

private:
        std::vector<Point> myPoints;
};
#endif
