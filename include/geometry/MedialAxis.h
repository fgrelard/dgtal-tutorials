#ifndef MEDIAL_AXIS_H
#define MEDIAL_AXIS_H

#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "geometry/DistanceToPointFunctor.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include <vector>

using namespace DGtal;

template <typename ImageFct, typename Point>
void checkPointForMedialAxis(const ImageFct& imageFct, std::vector<Point>& vPoints, const Point& p);

template <typename ImageFct, typename Point>
void checkPointForMedialAxis(const ImageFct& imageFct, std::vector<Point>& vPoints, const Point& p) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space,2> Distance;
	typedef DistanceToPointFunctor<Distance>         DistanceToPoint;
	typedef MetricAdjacency<Z3i::Space, 3>                Graph;
	typedef DistanceBreadthFirstVisitor< Graph, DistanceToPoint, std::set<Point> >
		DistanceVisitor;
	typedef typename DistanceVisitor::Node MyNode;
   
	double d = imageFct(p);
	Graph             graph;
	DistanceToPoint   d2pfct( Distance(), p );
	DistanceVisitor   visitor( graph, d2pfct, p );
	MyNode node;
	bool add = true;
	visitor.expand(); //Go to the first neighbour
	
	while ( ! visitor.finished() )
	{
		node = visitor.current(); // all the vertices of the same layer have been processed.
		
		if (node.second > 1) break; // we want to analyse only the direct neighbourhood (4- or 6- connexity)
		if (imageFct.domain().isInside(node.first) && node.second == 1) { //is inside domain
			float distanceToBoundary = imageFct(node.first);
			float minEnclosingRadius = sqrt(1+pow(d, 2));
			if (d <= 1 || distanceToBoundary >= minEnclosingRadius) {
				add = false;
			}
			visitor.expand();
		}
		else
			visitor.ignore();
	}
	if (add) {
		vPoints.push_back(p);
	}
}

#endif
 

