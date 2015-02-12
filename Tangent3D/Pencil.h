#ifndef PENCIL_H
#define PENCIL_H

#include <vector>
template <typename Point, typename Tangent, typename Vector>
class Pencil {
public:
	typedef Tangent T;
	typedef Vector Vector3d;
	typedef Point P;
public:
	Pencil(Point _point) : point(_point) {}
	Vector getTangent() const { return tangent; }
	Point getPoint() const { return point; }
	bool isUndefined() const { return undefined; }
	void compute_tangent(std::vector<Tangent> tangents);
private:
    Vector tangent;
	Point point;
	bool undefined = false;
};

#include "Pencil.ih"

#endif
