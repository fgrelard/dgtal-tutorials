#ifndef CURVE_ANALYZER_H
#define CURVE_ANALYZER_H

#include <vector>
#include <queue>

#include "graph/GraphEdge.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/MetricAdjacency.h"

namespace CurveAnalyzer {
	bool sameSet(const DGtal::Z3i::DigitalSet& first, const DGtal::Z3i::DigitalSet& second);

	double computeAngle(const DGtal::Z3i::Point& b,
						const DGtal::Z3i::Point& m,
						const DGtal::Z3i::Point& e);

	std::vector<GraphEdge*> neighboringEdges(const std::vector<GraphEdge*>& edges,
											 const DGtal::Z3i::DigitalSet& currentEdge,
											 const DGtal::Z3i::DigitalSet& branchingPoints);

	std::vector<GraphEdge*> hierarchicalDecomposition(const std::vector<DGtal::Z3i::DigitalSet>& edges,
													  const std::vector<DGtal::Z3i::Point>& endPoints,
													  const DGtal::Z3i::DigitalSet& branchingPoints);

	DGtal::Z3i::DigitalSet ensureConnexity(const DGtal::Z3i::DigitalSet& set);

	std::vector<DGtal::Z3i::Point> findEndPoints(const DGtal::Z3i::DigitalSet& set);


	template <typename Segmentation, typename Domain>
	DGtal::Z3i::DigitalSet dominantPointDetection(const Segmentation& segmentation,
												  const std::vector<DGtal::Z3i::Point>& skeletonOrdered,
												  const Domain& domain);


	DGtal::Z3i::DigitalSet detectCriticalPoints(const DGtal::Z3i::DigitalSet& skeleton);

	template <typename Segmentation, typename Domain>
	DGtal::Z3i::DigitalSet branchingPointDetection(const Segmentation& segmentation,
												   const std::vector<DGtal::Z3i::Point>& skeletonOrdered,
												   const Domain& domain);

	std::vector<DGtal::Z3i::DigitalSet> constructGraph(const std::vector<DGtal::Z3i::Point>& orderedCurve,
													   const DGtal::Z3i::DigitalSet& constraint);


	std::vector<DGtal::Z3i::Point> convertToOrientedEdge(const DGtal::Z3i::DigitalSet& edge);

	std::vector<DGtal::Z3i::Point> convertToOrientedEdge(const DGtal::Z3i::DigitalSet& edge, const DGtal::Z3i::Point& startingPoint);
};


bool CurveAnalyzer::sameSet(const DGtal::Z3i::DigitalSet& first, const DGtal::Z3i::DigitalSet& second) {
	if (first.size() != second.size()) return false;
	unsigned int cpt = 0;
	for (const DGtal::Z3i::Point& f : first) {
		if (second.find(f) != second.end())
			cpt++;
	}
	return (cpt == first.size());
}

double CurveAnalyzer::computeAngle(const DGtal::Z3i::Point& b,
								   const DGtal::Z3i::Point& m,
								   const DGtal::Z3i::Point& e) {
	DGtal::Z3i::RealVector v1 = (b - m).getNormalized();
	DGtal::Z3i::RealVector v2 = (e - m).getNormalized();
    double dot = v1.dot(v2);
	double angle = acos(dot);
	return angle;
}

std::vector<GraphEdge*> CurveAnalyzer::neighboringEdges(const std::vector<GraphEdge*>& edges,
									const DGtal::Z3i::DigitalSet& currentEdge,
									const DGtal::Z3i::DigitalSet& branchingPoints) {
	typedef DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> MetricAdjacency;

	std::vector<GraphEdge*> neighbors;
	DGtal::Z3i::Point branchPoint;
	for (const DGtal::Z3i::Point& b : branchingPoints) {
		if (currentEdge.find(b) != currentEdge.end() ) {
			branchPoint = b;
		}
	}

	std::vector<DGtal::Z3i::Point> nb;
	std::back_insert_iterator<std::vector<DGtal::Z3i::Point>> inserter(nb);
	MetricAdjacency::writeNeighbors(inserter, branchPoint);

	for (GraphEdge* edge : edges) {
		DGtal::Z3i::DigitalSet setEdge = edge->pointSet();
		if (sameSet(setEdge, currentEdge)) continue;
		for (const DGtal::Z3i::Point& n : nb) {
			if (setEdge.find(n) != setEdge.end())
				neighbors.push_back(edge);
		}
	}

	return neighbors;
}

std::vector<GraphEdge*> CurveAnalyzer::hierarchicalDecomposition(const std::vector<DGtal::Z3i::DigitalSet>& edges,
											const std::vector<DGtal::Z3i::Point>& endPoints,
											const DGtal::Z3i::DigitalSet& branchingPoints) {

	std::queue<GraphEdge*> edgeQueue;
	std::vector<GraphEdge*> hierarchyGraph;
	for (const DGtal::Z3i::DigitalSet& edge : edges) {
		GraphEdge* levelEdge = new GraphEdge(edge, std::numeric_limits<int>::max());

		for (const DGtal::Z3i::Point& e : endPoints) {
			if (edge.find(e) != edge.end()) {
				levelEdge->setLabel(1);
				edgeQueue.push(levelEdge);
				break;
			}
		}
		hierarchyGraph.push_back(levelEdge);
	}

	while (!edgeQueue.empty()) {
		GraphEdge* edgeCurrent  = edgeQueue.front();
		DGtal::Z3i::DigitalSet edgeCurrentSet = edgeCurrent->pointSet();
		edgeQueue.pop();
		std::vector<GraphEdge*> neighborEdges = neighboringEdges(hierarchyGraph, edgeCurrentSet, branchingPoints);
		for (GraphEdge* neighCurrent : neighborEdges) {
			int label = edgeCurrent->getLabel()+1;
			if (neighCurrent->getLabel() > label) {
				neighCurrent->setLabel(label);
				edgeQueue.push(neighCurrent);
			}
		}
	}

	return hierarchyGraph;

}

DGtal::Z3i::DigitalSet CurveAnalyzer::ensureConnexity(const DGtal::Z3i::DigitalSet& set) {
	DGtal::Z3i::DigitalSet cleanSet(set.domain());
	DGtal::Z3i::Object26_6 obj(DGtal::Z3i::dt26_6, set);
	DGtal::Z3i::DigitalSet & S = obj.pointSet();
	cleanSet = S;
	for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
		DGtal::Z3i::Object26_6 obj(DGtal::Z3i::dt26_6, cleanSet);
		if (obj.isSimple(*it)) {
		    cleanSet.erase(*it);
		}
	}

	return cleanSet;
}

std::vector<DGtal::Z3i::Point> CurveAnalyzer::findEndPoints(const DGtal::Z3i::DigitalSet& set) {
	DGtal::Z3i::Object26_6 objectSet(DGtal::Z3i::dt26_6, set);
	std::vector<DGtal::Z3i::Point> endPoints;
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		DGtal::Z3i::Point p = *it;
		std::vector<DGtal::Z3i::Point> neighbors;
		std::back_insert_iterator<std::vector<DGtal::Z3i::Point>> inserter(neighbors);
		objectSet.writeNeighbors(inserter, p);
		if (neighbors.size() <= 1)
			endPoints.push_back(p);

		//Is it in same quadrant: case connectivity != 26
		else {
			DGtal::Z3i::RealVector previous;
			bool isEndPoint = true;
			std::vector<DGtal::Z3i::RealVector> vectors;
			for (const DGtal::Z3i::Point& n : neighbors) {
				DGtal::Z3i::RealVector dir = (n - p).getNormalized();
				vectors.push_back(dir);
			}
			//Min angle (resp max dot product) determined by two points with one varying coordinate
			for (int i = 0; i < vectors.size(); i++) {
				for (int j = i+1; j < vectors.size(); j++) {
					if (vectors[i].dot(vectors[j]) <= (1/(1+sqrt(2))) )
						isEndPoint = false;
				}
			}
			if (isEndPoint)
				endPoints.push_back(p);
		}
	}
	return endPoints;
}



template <typename Segmentation, typename Domain>
DGtal::Z3i::DigitalSet CurveAnalyzer::dominantPointDetection(const Segmentation& segmentation,
									   const std::vector<DGtal::Z3i::Point>& skeletonOrdered,
									   const Domain& domain) {

	typedef typename Segmentation::SegmentComputerIterator DSS;

	DGtal::Z3i::DigitalSet dominantPoints(domain);
	DSS q = (segmentation.begin());
	DSS p = (++segmentation.begin());
	while (p != segmentation.end()) {
		DGtal::Z3i::Point endQ = *(--(q->end()));
		DGtal::Z3i::Point beginP = *(p->begin());
		int posP = find(skeletonOrdered.begin(), skeletonOrdered.end(), beginP) - skeletonOrdered.begin();
		int posQ = find(skeletonOrdered.begin(), skeletonOrdered.end(), endQ) - skeletonOrdered.begin();
	    if (posP < posQ) {
			DGtal::Z3i::Point beginQ = *(q->begin());
			DGtal::Z3i::Point endP = *(--(p->end()));
			double min = std::numeric_limits<double>::max();
			DGtal::Z3i::Point dominantCandidate;
			for (int i = posP; i < posQ; i++) {
				DGtal::Z3i::Point pointI = skeletonOrdered[i];
				double angle = computeAngle(beginQ, pointI, endP);
				if (angle < min) {
					min = angle;
					dominantCandidate = pointI;
				}
			}
			dominantPoints.insert(dominantCandidate);
		}
		++p;
		++q;
	}
	return dominantPoints;
}

template <typename Segmentation, typename Domain>
DGtal::Z3i::DigitalSet CurveAnalyzer::branchingPointDetection(const Segmentation& segmentation,
										const std::vector<DGtal::Z3i::Point>& skeletonOrdered,
										const Domain& domain) {

	typedef typename Segmentation::SegmentComputerIterator DSSIterator;
	typedef typename Segmentation::SegmentComputer DSS;

	DGtal::Z3i::DigitalSet branchingPoints(domain);

	for (const DGtal::Z3i::Point& s : skeletonOrdered) {

		std::vector<DSS> segments;

		for (DSSIterator p = segmentation.begin(), q = segmentation.end(); p != q; ++p) {
			DSS segment(*p);
			for (auto it = segment.begin(), ite = segment.end(); it != ite; ++it) {
				if (*it == s) {
					segments.push_back(segment);
					break;
				}
			}
		}

		int nb = 0;
		for (int i = 0; i < segments.size(); i++) {
			for (int j = i+1; j < segments.size(); j++) {
				DSS si = segments[i];
				DSS sj = segments[j];
				DGtal::Z3i::Point beginI = *(si.begin());
				DGtal::Z3i::Point endI = *(--(si.end()));
				DGtal::Z3i::Point beginJ = *(sj.begin());
				DGtal::Z3i::Point endJ = *(--(sj.end()));
				int posBI = find(skeletonOrdered.begin(), skeletonOrdered.end(), beginI) - skeletonOrdered.begin();
				int posEI = find(skeletonOrdered.begin(), skeletonOrdered.end(), endI) - skeletonOrdered.begin();
				int posBJ = find(skeletonOrdered.begin(), skeletonOrdered.end(), beginJ) - skeletonOrdered.begin();
				int posEJ = find(skeletonOrdered.begin(), skeletonOrdered.end(), endJ) - skeletonOrdered.begin();
				if (posBJ < posEI || posBI < posEJ) {
					nb++;
				}
			}
		}
		if (nb >= 2) {
			branchingPoints.insert(s);
		}
	}
	return branchingPoints;
}


std::vector<DGtal::Z3i::DigitalSet> CurveAnalyzer::constructGraph(const std::vector<DGtal::Z3i::Point>& orderedCurve,
									   const DGtal::Z3i::DigitalSet& constraint) {
	std::vector<DGtal::Z3i::DigitalSet> graph;
	int index = 0;
	graph.push_back(DGtal::Z3i::DigitalSet(constraint.domain()));
	DGtal::Z3i::Point previous;
	for (int i = 0, end = orderedCurve.size(); i < end; i++) {
		DGtal::Z3i::Point current = orderedCurve[i];
		if (DGtal::Z3i::l2Metric(previous,current) <= sqrt(3) || previous == DGtal::Z3i::Point())
			graph[index].insert(current);
		if (constraint.find(current) != constraint.end() || DGtal::Z3i::l2Metric(previous,current) > sqrt(3)) {
			index++;
			graph.push_back(DGtal::Z3i::DigitalSet(constraint.domain()));
			graph[index].insert(current);
		}
		previous = current;
	}
	return graph;
}

DGtal::Z3i::DigitalSet CurveAnalyzer::detectCriticalPoints(const DGtal::Z3i::DigitalSet& skeleton) {
	DGtal::Z3i::Object26_6 obj(DGtal::Z3i::dt26_6, skeleton);
	DGtal::Z3i::DigitalSet criticalPoints(skeleton.domain());
	for (const DGtal::Z3i::Point& s : skeleton) {
		std::vector<DGtal::Z3i::Point> neighbors;
		std::back_insert_iterator<std::vector<DGtal::Z3i::Point>> inserter(neighbors);
		obj.writeNeighbors(inserter, s);
		if (neighbors.size() > 2) {
			criticalPoints.insert(s);
		}
	}
	return criticalPoints;
}


std::vector<DGtal::Z3i::Point> CurveAnalyzer::convertToOrientedEdge(const DGtal::Z3i::DigitalSet& edge) {
	std::vector<DGtal::Z3i::Point> orientedEdge;
	if (edge.size() == 0) return orientedEdge;

	//DGtal::Z3i::DigitalSet thinEdge = ensureConnexity(edge);
	std::vector<DGtal::Z3i::Point> endPoints = findEndPoints(edge);
	DGtal::Z3i::Point start = *(endPoints.begin());
	return convertToOrientedEdge(edge, start);

}


std::vector<DGtal::Z3i::Point> CurveAnalyzer::convertToOrientedEdge(const DGtal::Z3i::DigitalSet& edge, const DGtal::Z3i::Point& startingPoint) {
	std::vector<DGtal::Z3i::Point> orientedEdge;
	if (edge.size() == 0) return orientedEdge;

	DGtal::Z3i::Object26_6 objEdge(DGtal::Z3i::dt26_6, edge);
	DGtal::Z3i::Point start = startingPoint;
	orientedEdge.push_back(start);
	bool toAdd = true;
	while (toAdd) {
		std::vector<DGtal::Z3i::Point> neighbors;
		std::back_insert_iterator<std::vector<DGtal::Z3i::Point>> inserter(neighbors);
		objEdge.writeNeighbors(inserter, start);
		unsigned int cpt = 0;
		for (const DGtal::Z3i::Point& n : neighbors) {
			if (std::find(orientedEdge.begin(), orientedEdge.end(), n) == orientedEdge.end()) {
				orientedEdge.push_back(n);
				start = n;
			}
			else
				cpt++;
		}
		if (cpt == neighbors.size())
			toAdd = false;
	}
	return orientedEdge;
}


#endif
