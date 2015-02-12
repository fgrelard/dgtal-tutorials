#ifndef MSTTANGENT_H
#define MSTTANGENT_H

template <typename Vector3D>
class MSTTangent {
public:
	typedef Vector3D Vector;
	MSTTangent() : position(0), v(Vector()) {}
	MSTTangent(int _position, Vector _v) : position(_position), v(_v) {}
	void computeEccentricity(float size) { eccentricity = position / size; }
private:
	int position;
public:
	Vector v;
	float eccentricity;
};

#endif
