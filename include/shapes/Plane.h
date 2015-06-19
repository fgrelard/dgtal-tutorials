#ifndef PLANE_H
#define PLANE_H

template <typename Point>
class Plane {
private:
	int myA;
	int myB;
	int myC;
	Point myOrigin;
public:
	Plane() : Plane(0, 0, 0, {0,0,0}) {}
	Plane(int a, int b, int c, Point origin) : myA(a), myB(b), myC(c), myOrigin(origin) {}
	Point getNormal() { return Point(myA, myB, myC); }
	Point getOrigin() { return myOrigin;}
};

#endif
