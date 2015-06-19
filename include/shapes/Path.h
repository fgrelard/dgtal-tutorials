#ifndef PATH_H
#define PATH_H

#include <vector>

template <typename Vertex>
class Path {
public:
	typedef typename std::vector<Vertex>::const_iterator Iterator;
public:
	Path() : myPath{0}, myIntersection{Vertex{}} {}
	Path(const std::vector<Vertex>& path, const Vertex& bel, const Vertex& intersection, const std::pair<Vertex, Vertex>& nearBel) : myPath{path}, myBel{bel}, myIntersection{intersection}, myNearBel{nearBel} {}
	Iterator begin() const { return myPath.cbegin(); }
	Iterator end() const { return myPath.cend(); }
	int size() const  { return myPath.size(); }
public:
	std::vector<Vertex> myPath;
	Vertex myBel;
	Vertex myIntersection;
	std::pair<Vertex, Vertex> myNearBel;

};
#endif
