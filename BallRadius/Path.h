#ifndef PATH_H
#define PATH_H

#include <vector>

template <typename Vertex>
class Path {
public:
	typedef typename std::vector<Vertex>::const_iterator Iterator;
public:
	Path(const std::vector<Vertex>& path, const Vertex& intersection) : myPath{path}, myIntersection{intersection} {}
	Iterator begin() const { return myPath.cbegin(); }
	Iterator end() const { return myPath.cend(); }
private:
	std::vector<Vertex> myPath;
	Vertex myIntersection;
};
#endif
