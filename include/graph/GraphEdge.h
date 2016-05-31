#ifndef GRAPH_EDGE_H
#define GRAPH_EDGE_H

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

class GraphEdge {
public:
	GraphEdge(DGtal::Z3i::DigitalSet points, int label) : myPoints(points), myLabel(label) {}
	int getLabel() const { return myLabel; }
	void setLabel(int label) { myLabel = label; }
	DGtal::Z3i::DigitalSet pointSet() const { return myPoints; }
private:
	DGtal::Z3i::DigitalSet myPoints;
	int myLabel;
};

#endif
