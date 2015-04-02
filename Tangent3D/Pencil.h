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
	Pencil(Point _point, Vector _tangent) : point(_point), tangent(_tangent) {}
	Vector getTangent() const { return tangent; }
	Point getPoint() const { return point; }
	bool isUndefined() const { return undefined; }
	void compute_tangent(std::vector<Tangent> tangents);
	friend inline bool operator<(const Pencil& it, const Pencil& other) {
		return it.point < other.point;
	}
private:
    Vector tangent;
	Point point;
	bool undefined = false;
};

#include "Pencil.ih"

#endif
