#ifndef SURFACE_TRAVERSAL_H
#define SURFACE_TRAVERSAL_H

#include <vector>
#include <map>
#include <set>
#include "DGtal/base/Common.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/helpers/StdDefs.h"
#include "geometry/DistanceToPointFunctor.h"

namespace SurfaceTraversal {
	/**
	 * A-star algorithm
	 */
	template <typename Point, typename Graph>
	std::vector<Point> AStarAlgorithm(const Graph& graph, const Point& source, const Point& destination);


	template <typename Point, typename Graph, typename Set>
	std::vector<Point> AStarAlgorithm(const Graph& graph, const Point& source, const Set& set);
	
	/**
	 * Dijkstra algorithm
	 */
	template <typename Point, typename Graph>
	std::vector<Point> DijkstraAlgorithm(const Graph& graph, const Point& source, const Point& destination);
	
	/**
	 * This function can be used for A star or Dijkstra depending on the visitor
	 */
	template <typename Point, typename Visitor>
	std::vector<Point> shortestPath(Visitor& visitor, const Point& source, const Point& goal);


	template <typename Point, typename Visitor, typename Set>
	std::vector<Point> shortestPath(Visitor& visitor, const Point& source, const Set& set);
	
	/**
	 * Reconstructs a path
	 */
	template <typename Point>
	std::vector<Point> reconstructPath(const std::map<Point, Point>& aMapPrevious, const Point& source, const Point& goal);
}

template <typename Point, typename Graph>
std::vector<Point> SurfaceTraversal::AStarAlgorithm(const Graph& graph, const Point& source, const Point& destination) {
	typedef DistanceToPointFunctor<DGtal::Z3i::L2Metric> Distance;
	typedef DGtal::DistanceBreadthFirstVisitor<Graph, Distance, std::set<Point>> Visitor;
	Distance distance(DGtal::Z3i::L2Metric(), source);
	Visitor visitor(graph, distance, source);
	std::vector<Point> path = SurfaceTraversal::shortestPath(visitor, source, destination);
	return path;
}

template <typename Point, typename Graph, typename Set>
std::vector<Point> SurfaceTraversal::AStarAlgorithm(const Graph& graph, const Point& source, const Set& set) {
	typedef DistanceToPointFunctor<DGtal::Z3i::L2Metric> Distance;
	typedef DGtal::DistanceBreadthFirstVisitor<Graph, Distance, std::set<Point>> Visitor;
	Distance distance(DGtal::Z3i::L2Metric(), source);
	Visitor visitor(graph, distance, source);
	std::vector<Point> path = SurfaceTraversal::shortestPath(visitor, source, set);
	return path;
}

template <typename Point, typename Graph>
std::vector<Point> SurfaceTraversal::DijkstraAlgorithm(const Graph& graph, const Point& source, const Point& destination) {
	typedef DGtal::BreadthFirstVisitor<Graph, std::set<Point>> Visitor;
	Visitor visitor(graph, source);
	std::vector<Point> path = SurfaceTraversal::shortestPath(visitor, source, destination);
	return path;
}

template <typename Point, typename Visitor>
std::vector<Point> SurfaceTraversal::shortestPath(Visitor& visitor, const Point& source, const Point& goal) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	std::vector<Point> thePath;
	std::map<Point, Point> aMapPrevious;
	while (!visitor.finished()) {
		node = visitor.current();
		if (node.first == goal) {
			return reconstructPath(aMapPrevious, source, goal);
		}
		std::vector<Point> neighbors;
		std::back_insert_iterator<std::vector<Point>> iter(neighbors);
		visitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			auto itMapExisting = aMapPrevious.find(*it);
			if (itMapExisting == aMapPrevious.end()) {
				aMapPrevious[*it] = node.first;
			} 
		}
		visitor.expand();
	}
	return thePath; //path not found = empty vector
}


template <typename Point, typename Visitor, typename Set>
std::vector<Point> SurfaceTraversal::shortestPath(Visitor& visitor, const Point& source, const Set& set) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	std::vector<Point> thePath;
	std::map<Point, Point> aMapPrevious;
	while (!visitor.finished()) {
		node = visitor.current();
		if (source != node.first && set.find(node.first) != set.end()) {
			return reconstructPath(aMapPrevious, source, node.first);
		}
		std::vector<Point> neighbors;
		std::back_insert_iterator<std::vector<Point>> iter(neighbors);
		visitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			auto itMapExisting = aMapPrevious.find(*it);
			if (itMapExisting == aMapPrevious.end()) {
				aMapPrevious[*it] = node.first;
			}
		}
		visitor.expand();
	}
	return thePath; //path not found = empty vector
}


template <typename Point>
std::vector<Point> SurfaceTraversal::reconstructPath(const std::map<Point, Point>& aMapPrevious, const Point& source, const Point& goal) {
	std::vector<Point> path;
	Point aPoint = goal;
	
	while (aPoint != source) {
		path.push_back(aPoint);
		aPoint = aMapPrevious.at(aPoint);
	}
	
	
	path.push_back(source);
	return path;
}

#endif
