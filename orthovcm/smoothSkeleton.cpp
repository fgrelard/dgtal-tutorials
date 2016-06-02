/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file dvcm-2d.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/04/15
 *
 * Computes the Voronoi Covariance Measure of a list of 2D digital
 * points. Displays the resulting normal vector and feature detection.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <iterator>
#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "geometry/VoronoiCovarianceMeasure.h"
#include "geometry/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/math/Histogram.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/graph/GraphVisitorRange.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/TangentUtils.h"
#include "geometry/SliceUtils.h"
#include "geometry/Pencil.h"
#include "geometry/MSTTangent.h"
#include "geometry/PointUtil.h"
#include "geometry/WeightedPoint.h"
#include "geometry/MedialAxis.h"
#include "geometry/ImageUtil.h"
#include "surface/SurfaceUtils.h"
#include "surface/Morphomaths.h"
#include "clustering/diana.hpp"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "geometry/WeightedPointCount.h"
#include "geometry/SaddleComputer.h"
#include "surface/SurfaceTraversal.h"
#include "Statistics.h"
#include "shapes/Ball.h"
///////////////////////////////////////////////////////////////////////////////
#include "ReverseHatPointFunction.h"


#include "graph/GraphEdge.h"
#include "geometry/CurveAnalyzer.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;



Z3i::Point extractNearestNeighborInSetFromPoint(const Z3i::DigitalSet& aSet, const Z3i::RealPoint& aPoint) {
	double distanceMin = numeric_limits<double>::max();
	Z3i::Point toReturn;
	for (auto it = aSet.begin(), ite = aSet.end(); it != ite; ++it) {
		double distanceToPoint = sqrt(pow(((*it)[0] - aPoint[0]), 2) + pow(((*it)[1] - aPoint[1]), 2) + pow(((*it)[2] - aPoint[2]), 2));
		if (distanceToPoint < distanceMin || (distanceMin == distanceToPoint && aPoint > *it)) {
			distanceMin = distanceToPoint;
			toReturn = *it;
		}
	}
	return toReturn;
}

Z3i::DigitalSet reduceClustersToCenters(const Z3i::DigitalSet& clusters) {
	Z3i::DigitalSet newClusters(clusters.domain());

	Z3i::Object26_6 obj(Z3i::dt26_6, clusters);
	vector<Z3i::Object26_6> objects;
	back_insert_iterator<vector<Z3i::Object26_6>> inserter(objects);
    obj.writeComponents(inserter);
	for (const Z3i::Object26_6& o : objects) {
		Z3i::DigitalSet currentPointSet = o.pointSet();
		if (currentPointSet.size() > 1) {
			Z3i::RealPoint center = Statistics::extractCenterOfMass3D(currentPointSet);
			Z3i::Point centerInSet = extractNearestNeighborInSetFromPoint(currentPointSet, center);
			newClusters.insert(centerInSet);
		}
		else {
			newClusters.insert(currentPointSet.begin(), currentPointSet.end());
		}
	}
	return newClusters;
}

template <typename DTL2>
Z3i::Point findMaxDTInSet(const Z3i::DigitalSet& set, const DTL2 dt, const Z3i::Point& junctionPoint) {
	double maxDT = 0.0;
	Z3i::Point maxDTPoint;
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		Z3i::Point point = *it;
		if (dt(point) > maxDT) {
			maxDT = dt(point);
			maxDTPoint = point;
		}
		else if (dt(point) == maxDT) {
			if (Z3i::l2Metric(junctionPoint, point) < Z3i::l2Metric(junctionPoint, maxDTPoint))
				maxDTPoint = point;
		}
	}
	return maxDTPoint;
}

template <typename Image>
double computeRadiusFromIntersection(const Image& volume, const Z3i::Point& point, const Z3i::RealPoint& normal,
									 double radius, int column, int scalar=1) {
	typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, typename Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;
	typedef ImageSelector<Z2i::Domain, bool>::Type Image2D;
	DGtal::functors::Identity idV;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-radius, -radius, -radius), volume.domain().upperBound() + Z3i::Point(radius, radius, radius));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0),
									  DGtal::Z2i::Point(radius, radius));


	DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, point, normal, radius, domain3Dyup.lowerBound());

	ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
	Image2D processImage = ImageUtil::convertImage<Image2D>(extractedImage);
	Z2i::Point center(radius/2, radius/2);
	Z2i::DigitalSet aSet = VCMUtil::extractConnectedComponent(processImage, center, 1, 255);
	Eigen::MatrixXd covmatrix = Statistics::computeCovarianceMatrix<Eigen::MatrixXd>(aSet);
	if (covmatrix.size() == 0) return 0;
	Z2i::RealVector projection = Statistics::extractEigenVector<Z2i::RealVector>(covmatrix, column);
	Z2i::Point trackedPoint = PointUtil::trackPoint(center, aSet, projection*scalar);
	double distance = Z2i::l2Metric(center, trackedPoint);
	return distance;
}



template <typename DTL2>
Z3i::DigitalSet extractSurfacePoints(const Z3i::DigitalSet& intersection,
									 const DTL2& dt) {
	Z3i::DigitalSet surfacePoints(intersection.domain());
	for (const Z3i::Point& p : intersection) {
		if (dt(p) <= 1)
			surfacePoints.insert(p);

	}
	return surfacePoints;
}

Z3i::DigitalSet extractSetFromRadiusAndCenter(const Z3i::DigitalSet& intersection,
											  const Z3i::Point& center,
											  double radius) {
	Ball<Z3i::Point> ball(center, radius);
	vector<Z3i::Point> restrictedIntersection = ball.intersection(intersection);
	Z3i::DigitalSet restrictedSet(intersection.domain());
	restrictedSet.insert(restrictedIntersection.begin(), restrictedIntersection.end());
	return restrictedSet;
}

template <typename DTL2>
Z3i::DigitalSet extractRelevantSurfacePointsForVolVCM(const Z3i::DigitalSet& intersection,
									 const DTL2& dt,
									 const Z3i::Point& center) {
	typedef Z3i::Object26_6 ObjectType;

	double maxSize = intersection.size();
	double currentSize = 0;
	unsigned int nbConnectedComponents = 0;
	double radius = 1.;
	Z3i::DigitalSet surfels(intersection.domain());

	while (currentSize < maxSize && nbConnectedComponents != 1) {
		Z3i::DigitalSet plane = extractSetFromRadiusAndCenter(intersection, center, radius);
	    surfels = extractSurfacePoints(intersection, dt);
		ObjectType objectImage(Z3i::dt26_6, surfels);
		vector<ObjectType> objects;
		back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	    nbConnectedComponents = objectImage.writeComponents(inserter);
		currentSize = plane.size();
		radius++;
	}
	return surfels;
}

template <typename DTL2>
Z3i::Point projectionPoint(const DTL2& dt,
						   const Z3i::Point& center,
						   const Z3i::RealVector& direction) {
	Z3i::Point proj = center;
	int scalar = 1;
	while (dt.domain().isInside(proj) && dt(proj) > 1) {
		scalar++;
		proj = center + direction*scalar;
	}
	return proj;
}

template <typename DTL2>
Z3i::RealVector signVector(const DTL2& dt,
						   const Z3i::Point& center,
						   const Z3i::RealVector& direction) {
	Z3i::Point projCenter = projectionPoint(dt, center, direction);
	Z3i::Point otherProjCenter = projectionPoint(dt,center, -direction);

	if (Z3i::l2Metric(projCenter, center) > Z3i::l2Metric(otherProjCenter, center))
		return -direction;
	return direction;
}

template <typename DTL2>
Z3i::DigitalSet extractSurfacePointsWithEigenVectors(const Z3i::DigitalSet& intersection, const DTL2& dt,
													 const Z3i::Point& center, double radius) {
	typedef Eigen::MatrixXd MatrixXd;

	Z3i::DigitalSet surfelSet(intersection.domain());
	MatrixXd covmat = Statistics::computeCovarianceMatrix<MatrixXd>(intersection);
	if (covmat.size() == 0) return surfelSet;
	Z3i::RealVector v2 = Statistics::extractEigenVector<Z3i::RealVector>(covmat, 2);
	v2 = signVector(dt, center, v2);
	v2 = -v2;
	Z3i::Point projectionCenter = center;
	int scalar=1;
	while (dt.domain().isInside(projectionCenter) && dt(projectionCenter) >= 1) {
		projectionCenter = center - v2*scalar;
		scalar++;
	}
	Ball<Z3i::Point> ball(projectionCenter, radius);
	vector<Z3i::Point> surfels = ball.intersection(intersection);
	surfelSet.insert(surfels.begin(), surfels.end());
	return surfelSet;
}



template <typename VCM, typename KernelFunction, typename Container>
double distanceToDelineateSubVolume(const Z3i::Point& current,
									const Z3i::Point& b,
									const vector<GraphEdge*> neighborEdge,
									const Container& setVolume,
									const Z3i::DigitalSet& saddlePoints) {

	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;
	double distance = 0;
	for (GraphEdge* edge : neighborEdge) {
		Z3i::DigitalSet setEdge = edge->pointSet();
		if (setEdge.find(current) != setEdge.end()) continue;
		for (const Z3i::Point& e : setEdge) {
			double radius = setEdge.size() * 0.4;
			VCM vcm(20, ceil(radius), l2, false);
			vcm.init(setEdge.begin(), setEdge.end());
			KernelFunction chi(1.0, radius);
			Z3i::RealPoint normal;
			Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, saddlePoints, setVolume,
																				 e, normal,
																				 0, radius, radius*2, false);
			for (const Z3i::Point& s : saddlePoints) {
				if (connectedComponent3D.find(s) != connectedComponent3D.end()) {
					double currentDistance = Z3i::l2Metric(e, b);
					if (currentDistance > distance) {
						distance = currentDistance;
					}
				}

			}
		}
	}
	return distance;
}

template <typename VCM, typename KernelFunction>
map<Z3i::Point, Z3i::RealPoint> computePlanesForSubVolume(const Z3i::Point& b,
														  const vector<GraphEdge*> neighborEdge,
														  double distance) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;

	map<Z3i::Point, Z3i::RealPoint> normals;
	for (GraphEdge* edge: neighborEdge) {
		Z3i::DigitalSet setEdge = edge->pointSet();
		double maxDifference = std::numeric_limits<double>::max();
		Z3i::Point candidate;
		for (const Z3i::Point& e : setEdge) {
			double difference = std::abs(Z3i::l2Metric(e, b) - distance);
			if (difference < maxDifference) {
				candidate  = e;
				maxDifference = difference;
			}
		}

		double radius = setEdge.size() * 0.4;
		VCM vcm(20, ceil(radius), l2, false);
		vcm.init(setEdge.begin(), setEdge.end());
		KernelFunction chi(1.0, radius);
		Z3i::RealPoint normal = VCMUtil::computeNormalFromVCM(candidate, vcm, chi, 0);
		Z3i::RealPoint directionVector = (candidate - b).getNormalized();
		normal = (normal.dot(directionVector) < 0) ? normal : -normal;
		normals[candidate] = normal;
	}
	return normals;
}

template <typename VCM, typename KernelFunction, typename Container>
Container constructSubVolume(const Z3i::Point& current,
							 const vector<GraphEdge*>& graph,
							 const Z3i::DigitalSet& currentEdge,
							 const Container& setVolume,
							 const Z3i::DigitalSet& branchingPoints,
							 const Z3i::DigitalSet& saddlePoints) {
	typedef MetricAdjacency<Z3i::Space, 3> MAdj;

	vector<Z3i::Point> neighbors;
	back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
	MAdj::writeNeighbors(inserter, current);

	Z3i::Point b;
	for (const Z3i::Point& n : neighbors) {
		if (branchingPoints.find(n) != branchingPoints.end())
			b = n;
	}
	if (b == Z3i::Point())
		return setVolume;

	vector<GraphEdge*> neighborEdge = CurveAnalyzer::neighboringEdges(graph,
													   currentEdge,
													   branchingPoints);

	double distance = distanceToDelineateSubVolume<VCM, KernelFunction>(current, b, neighborEdge, setVolume, saddlePoints);
	map<Z3i::Point, Z3i::RealPoint> normals = computePlanesForSubVolume<VCM, KernelFunction>(b, neighborEdge, distance);

	Container subVolume;
	for (auto it = setVolume.begin(), ite = setVolume.end();
		 it != ite; ++it) {
		Z3i::Point pointVolume = (*it)->myPoint;
		int cpt = 0;
		for (const auto& pair : normals) {
			if (pair.first == current)
				continue;
			if (VCMUtil::abovePlane(pointVolume, pair.second, pair.first)) {
				cpt++;
			}
		}
		bool toAdd = (cpt >= normals.size()-1);
		if (toAdd)
			subVolume.insert(*it);
	}
	return subVolume;
}



template <typename VCM, typename KernelFunction, typename Container>
Z3i::DigitalSet constructSubVolumeWithTangent(const Z3i::Point& current,
										const vector<GraphEdge*>& graph,
										const Z3i::DigitalSet& currentEdge,
										const Container& setVolume,
										const Z3i::DigitalSet& branchingPoints,
										const Z3i::DigitalSet& saddlePoints,
										const Z3i::DigitalSet& setSurface) {
	typedef MetricAdjacency<Z3i::Space, 3> MAdj;
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;

	Z3i::DigitalSet digitalSetVolume(saddlePoints.domain());
	for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
		digitalSetVolume.insert((*it)->myPoint);
	}

	vector<Z3i::Point> neighbors;
	back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
	MAdj::writeNeighbors(inserter, current);

	Z3i::Point b;
	for (const Z3i::Point& n : neighbors) {
		if (branchingPoints.find(n) != branchingPoints.end())
			b = n;
	}
	if (b == Z3i::Point())
		return digitalSetVolume;


	vector<GraphEdge*> neighborEdge = CurveAnalyzer::neighboringEdges(graph,
													   currentEdge,
													   branchingPoints);
	double radius = currentEdge.size() * 0.4;
	VCM vcm(20, ceil(radius), l2, false);
	vcm.init(currentEdge.begin(), currentEdge.end());
	KernelFunction chi(1.0, radius);
	Z3i::RealPoint normal = VCMUtil::computeNormalFromVCM(current, vcm, chi, 0);

	Z3i::DigitalSet subVolume(digitalSetVolume.domain());
	Z3i::Object26_6 obj(Z3i::dt26_6, setSurface);
	for (GraphEdge* edge : neighborEdge) {
		Z3i::DigitalSet setEdge = edge->pointSet();
		for (const Z3i::Point& e : setEdge) {
			Z3i::Point trackedPoint = PointUtil::trackPoint(e, digitalSetVolume, normal);
			Z3i::Point trackedPointMinus = PointUtil::trackPoint(e, digitalSetVolume, -normal);
			if (saddlePoints.find(trackedPoint) != saddlePoints.end() ||
				saddlePoints.find(trackedPointMinus) != saddlePoints.end()) {
				vector<Z3i::Point> path = SurfaceTraversal::AStarAlgorithm(obj, trackedPoint, trackedPointMinus);
				subVolume.insert(path.begin(), path.end());
				break;
			}
		}
	}
	return subVolume;
}


std::set<Z3i::Point> vectorsToCompare(const std::vector<Z3i::Point>& directionVectors) {

	set<Z3i::Point> vectorsToCompare;
	for (const Z3i::Point& directionVector: directionVectors) {
		for (const Z3i::Point& otherDirectionVector : directionVectors) {
			if (directionVector == otherDirectionVector) continue;
			if (directionVector.dot(otherDirectionVector) > 0) {
				vectorsToCompare.insert(directionVector);
				vectorsToCompare.insert(otherDirectionVector);
			}
		}
	}
	return vectorsToCompare;
}

std::set<Z3i::Point> vectorsToPoint(const std::set<Z3i::Point>& vectorsToCompare,
									const Z3i::Point& branchingPoint) {
	set<Z3i::Point> vectorsToPoint;
	for (const Z3i::Point& p : vectorsToCompare) {
		vectorsToPoint.insert(branchingPoint + p);
	}
	return vectorsToPoint;
}


vector<GraphEdge*> edgesAssociatedWithPoints(const std::set<Z3i::Point>& point,
											 const vector<GraphEdge*>& graph) {
	vector<GraphEdge*> edges;
	for (GraphEdge* edge : graph) {
		Z3i::DigitalSet setEdge = edge->pointSet();
		for (const Z3i::Point& p : point) {
			if (setEdge.find(p) != setEdge.end())
				edges.push_back(edge);
		}
	}
	return edges;
}

template <typename VCM, typename KernelFunction, typename Container, typename Domain>
map<Z3i::Point, vector<pair<Z3i::Point, Z3i::RealPoint> > > constructSubVolumeWithPlanes(const vector<GraphEdge*>& graph,
															 const Container& setVolume,
															 const Domain& domain,
															 const vector<Z3i::Object26_6>& saddlePoints) {
	typedef MetricAdjacency<Z3i::Space, 3> MAdj;
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;

	map<Z3i::Point, vector<pair<Z3i::Point, Z3i::RealPoint>>> saddlesToPlanes;
	for (const Z3i::Object26_6& o : saddlePoints) {
		Z3i::Point saddleCC = *(o.pointSet().begin());
		saddlesToPlanes[saddleCC] = vector<pair<Z3i::Point, Z3i::RealPoint>>();
	}

	Metric l2;
	int cpt = 0;

	for (GraphEdge* edge : graph) {
		Z3i::DigitalSet setEdge = edge->pointSet();
		if (setEdge.size() <= 2) continue;
		bool found = false;
		for (const Z3i::Point& p : setEdge) {
			if (found) break;
			double radius = setEdge.size() * 0.4;
			VCM vcm(20, ceil(radius), l2, false);
			vcm.init(setEdge.begin(), setEdge.end());
			KernelFunction chi(1.0, radius);
			Z3i::RealPoint normal;
			Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domain, setVolume,
																				 p, normal,
																				 0, radius, radius*2, false);
			for (const Z3i::Object26_6& o : saddlePoints) {
				Z3i::DigitalSet saddleCC = o.pointSet();
				if (found) break;
				for (const Z3i::Point& sp : saddleCC) {
					if (connectedComponent3D.find(sp) != connectedComponent3D.end()) {

						saddlesToPlanes[*saddleCC.begin()].push_back(make_pair(p, normal));
						found  = true;
						break;
					}
				}
			}
		}
	}

	for (auto it = saddlesToPlanes.begin(), ite = saddlesToPlanes.end(); it != ite; ++it) {
		vector<pair<Z3i::Point, Z3i::RealPoint>> pointToVectors = it->second;
		trace.info() << pointToVectors.size() << endl;
		if (pointToVectors.size() == 2) {
			Z3i::Point p1 = (*(pointToVectors.begin())).first;
			Z3i::RealPoint p2 =  (*(++pointToVectors.begin())).first;
			Z3i::RealPoint n1 = (*(pointToVectors.begin())).second;
			Z3i::RealPoint n2 =  (*(++pointToVectors.begin())).second;
			Z3i::RealPoint nrot = n1.crossProduct(n2);
			Z3i::RealPoint nrot2 = n2.crossProduct(n1);
			double alpha = acos(n1.dot(n2));
			float angleRotation = std::abs(M_PI/2 - alpha);
			Eigen::Affine3f t(Eigen::AngleAxisf(angleRotation, Eigen::Vector3f(nrot[0], nrot[1], nrot[2])));
			Eigen::Affine3f t2(Eigen::AngleAxisf(angleRotation, Eigen::Vector3f(nrot2[0], nrot2[1], nrot2[2])));
			Eigen::Vector3f n1e(n1[0], n1[1], n1[2]);
			Eigen::Vector3f n2e(n2[0], n2[1], n2[2]);
			Eigen::Vector3f rotn1e = t.linear() * n1e;
			Eigen::Vector3f rotn2e = t2.linear() * n2e;
			Eigen::Vector3f otherRotn1e = t2.linear() * n1e;
			Eigen::Vector3f otherRotn2e = t.linear() * n2e;
			Z3i::RealPoint rotn1(rotn1e[0], rotn1e[1], rotn1e[2]);
			Z3i::RealPoint rotn2(rotn2e[0], rotn2e[1], rotn2e[2]);
			Z3i::RealPoint otherRotn1(otherRotn1e[0], otherRotn1e[1], otherRotn1e[2]);
			Z3i::RealPoint otherRotn2(otherRotn2e[0], otherRotn2e[1], otherRotn2e[2]);
			if (std::abs(otherRotn1.dot(otherRotn2)) < std::abs(rotn1.dot(rotn2))) {
				rotn1 = otherRotn1;
				rotn2 = otherRotn2;
			}

			pointToVectors.clear();
			pointToVectors.push_back(make_pair(p1, rotn1));
			pointToVectors.push_back(make_pair(p2, rotn2));
			saddlesToPlanes[it->first] = pointToVectors;
		}
		else {
			saddlesToPlanes[it->first].clear();
		}
	}


	return saddlesToPlanes;
}

template <typename VCM, typename KernelFunction, typename Container>
vector<pair<Z3i::Point, double>> areaProfile(const Z3i::DigitalSet& setEdge, const Container& setVolume) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;
   
	auto domain = setEdge.domain();
	vector<Z3i::Point> orientedEdge = CurveAnalyzer::convertToOrientedEdge(setEdge);

	vector<pair<Z3i::Point, double>> areasEdge;
    auto itP = orientedEdge.begin();
	while (itP != orientedEdge.end()) {
		Z3i::Point p = *itP;
		double radius = setEdge.size() * 0.4;
		VCM vcm(20, ceil(radius), l2, false);
		vcm.init(setEdge.begin(), setEdge.end());
		KernelFunction chi(1.0, radius);
		Z3i::RealPoint normal;
		Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domain, setVolume,
																			 p, normal,
																			 0, radius, radius*2, false);

		double area = connectedComponent3D.size();
		areasEdge.push_back(make_pair(*itP, area));
		++itP;
	}
	return areasEdge;
}


template <typename VCM, typename KernelFunction, typename Domain, typename WeightedContainer>
pair<Z3i::Point, Z3i::RealPoint> pointToNormal(const Z3i::Point&  point, const Z3i::DigitalSet& setVCM, const Domain& domain,
											  const WeightedContainer& setVolume) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;

	map<Z3i::Point, Z3i::RealPoint> aMap;
	double radius = setVCM.size() * 0.4;
	radius = (radius < 2) ? 2 : radius;
	VCM vcm(20, ceil(radius), l2, false);
	vcm.init(setVCM.begin(), setVCM.end());
	KernelFunction chi(1.0, radius);
	Z3i::RealPoint normal;
	Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domain, setVolume,
																		 point, normal,
																		 0, radius, radius*2, false);
	return make_pair(point, normal);
	
}

template <typename VCM, typename KernelFunction, typename Domain, typename WeightedContainer>
Z3i::DigitalSet associatedPlane(const Z3i::Point& point, const Z3i::DigitalSet& setVCM, const Domain& domain,
											  const WeightedContainer& setVolume) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;

	map<Z3i::Point, Z3i::RealPoint> aMap;
	double radius = setVCM.size() * 0.4;
	VCM vcm(20, ceil(radius), l2, false);
	vcm.init(setVCM.begin(), setVCM.end());
	KernelFunction chi(1.0, radius);
	Z3i::RealPoint normal;
	Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domain, setVolume,
																		 point, normal,
																		 0, radius, radius*2, false);
	return connectedComponent3D;
}

vector<Z3i::Point> pointsAboveAverage(const map<Z3i::Point, double>& pointToAreas) {
	vector<double> areas;
	for (auto it = pointToAreas.begin(), ite = pointToAreas.end(); it != ite; ++it) {
		areas.push_back(it->second);
	}

	vector<Z3i::Point> points;
	double mean = Statistics::mean(areas);
	for (auto it = pointToAreas.begin(), ite = pointToAreas.end(); it!=ite; ++it) {
		Z3i::Point currentPoint = it->first;
		double currentArea = it->second;
		if (currentArea > mean) {
			points.push_back(currentPoint);
		}
	}
	return points;
}

Z3i::Point pointsVaryingNeighborhood(const vector<pair<Z3i::Point, double>>& pointToAreas) {
	typedef MetricAdjacency<Z3i::Space, 3> MAdj;
	
	Z3i::Point candidate;
	double previousFactor = 0;
	for (const auto& pair : pointToAreas) {
		Z3i::Point p = pair.first;
		double currentValue = pair.second;
		double radius = sqrt(currentValue / M_PI);
		double currentArea = currentValue + 2 * radius * M_PI + M_PI;
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
	    MAdj::writeNeighbors(inserter, p);
		for (const Z3i::Point& n : neighbors) {
			auto nInMap = find_if(pointToAreas.begin(),pointToAreas.end(), [&](const std::pair<Z3i::Point, double>& pair) {
					return pair.first == n;
				});
			if (nInMap != pointToAreas.end()) {
				double valueNeighbor = nInMap->second;
				double factor = valueNeighbor / currentArea;
				if (factor > previousFactor) {
					previousFactor = factor;
					candidate = p;
				}				
			}
		}	   
	}
	return candidate;
}


vector<Z3i::DigitalSet> computeSubVolumes(const Z3i::DigitalSet& volume,
								   const vector<Z3i::DigitalSet>& planes) {
	vector<Z3i::DigitalSet> subVolumes;

	subVolumes.push_back(volume);
	for (const Z3i::DigitalSet& plane : planes) {
		for (auto it = subVolumes.begin(), ite = subVolumes.end(); it!=ite; ++it) {
			Z3i::DigitalSet subVol = *it;
			set<Z3i::Point> subVolSet(subVol.begin(), subVol.end());
			set<Z3i::Point> planeSet(plane.begin(), plane.end());
			
			set<Z3i::Point> difference;
			set_difference(subVolSet.begin(), subVolSet.end(),
						   planeSet.begin(), planeSet.end(),
						   std::inserter(difference, difference.end()));
			Z3i::DigitalSet diffDigitalSet(subVol.domain());
			diffDigitalSet.insert(difference.begin(), difference.end());
			Z3i::Object6_26 objDiff(Z3i::dt6_26, diffDigitalSet);
			vector<Z3i::Object6_26> cc;
			back_insert_iterator< vector<Z3i::Object6_26> > iterator( cc );
			unsigned nbcc = objDiff.writeComponents(iterator);
			if (nbcc >= 2) {
				subVolumes.erase(it);
				for (auto itCC = cc.begin(), itCCe = cc.end();
					 itCC != itCCe; ++itCC) {
					Z3i::DigitalSet part = itCC->pointSet();
//					part.insert(plane.begin(), plane.end());
					subVolumes.push_back(part);
				}
				break;
			}
		}
	}
	return subVolumes;
}


vector<Z3i::DigitalSet> delineatingPlanes(const Z3i::DigitalSet& subVolume, const vector<Z3i::DigitalSet>& setOfAllPlanes) {
	vector<Z3i::DigitalSet> delineating;
	for (const Z3i::DigitalSet& plane : setOfAllPlanes) {
		Z3i::DigitalSet dilatedPlane(plane.domain());
		for (const Z3i::Point& pointPlane : plane) {
			dilatedPlane.insert(pointPlane);
			vector<Z3i::Point> neighbors;
			back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
			MetricAdjacency<Z3i::Space, 3>::writeNeighbors(inserter, pointPlane);
			for (const Z3i::Point& n : neighbors) {
				dilatedPlane.insert(n);
			}
		}
		for (const Z3i::Point& p : subVolume) {
				
			if (dilatedPlane.find(p) != dilatedPlane.end()) {
				delineating.push_back(plane);
				break;
			}
		}
	}
	return delineating;
}


vector<Z3i::DigitalSet> twoClosestPlanesHausdorff(const vector<Z3i::DigitalSet>& delineatingPlanes) {
	double distanceMin = numeric_limits<double>::max();
	vector<Z3i::DigitalSet> pairPlanes;
	for (size_t i = 0; i < delineatingPlanes.size(); i++) {
		for (size_t j = i+1; j < delineatingPlanes.size(); j++) {
			Z3i::DigitalSet first = delineatingPlanes[i];
			Z3i::DigitalSet second = delineatingPlanes[j];
			double dh = hausdorffDistance(first, second);
			if (dh < distanceMin) {
				vector<Z3i::DigitalSet> currentPair;
				currentPair.push_back(first);
				currentPair.push_back(second);
				pairPlanes = currentPair;
			}
		}
	}
	return pairPlanes;
}


template <typename VCM, typename KernelFunction, typename Container>
void areaProfileToStdout(const vector<Z3i::DigitalSet>& edgeGraph,
						 const Container& setVolumeWeighted) {
	int loop = 0;
	for (const Z3i::DigitalSet& edge : edgeGraph) {
		if (edge.size() < 10) continue;
		vector<pair<Z3i::Point, double>> pointToAreas = areaProfile<VCM, KernelFunction>(edge, setVolumeWeighted);
		int j = 0;
		for (const auto& pair : pointToAreas) {			
			cout << j << " " << pair.second << endl;
			j++;
		}
		loop++;
	}
}


///////////////////////////////////////////////////////////////////////////////
int main( int  argc, char**  argv )
{


	typedef Z3i::Space Space;
	typedef Z3i::Point Point;
	typedef Z3i::RealPoint RealPoint;
	typedef Z3i::RealVector RealVector;
	typedef HyperRectDomain<Space> Domain;
	typedef ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric
	typedef EigenDecomposition<3,double> LinearAlgebraTool;
	typedef LinearAlgebraTool::Matrix Matrix;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef VoronoiCovarianceMeasure<Space,Metric> VCM;


	typedef MSTTangent<Point> Tangent;
	typedef Pencil<Point, Tangent, RealPoint> Pencil;

	typedef WeightedPoint<Z3i::RealPoint> WeightedRealPoint;
	typedef WeightedPoint<Z3i::Point> WeightedPoint;
	typedef MetricAdjacency<Space, 3> MetricAdjacency;
	typedef WeightedPointCount<Point> WeightedPointCount;
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef DGtal::DistanceTransformation<Space, ThresholdedImage, Z3i::L2Metric> DTL2;
	typedef functors::BallConstantPointFunction<Point, double> KernelFunction;

	typedef Z3i::KSpace KSpace;
	typedef DGtal::functors::NotPointPredicate<ThresholdedImage> BackgroundPredicate;
	typedef ImplicitDigitalSurface< KSpace, BackgroundPredicate > DigitalSurfaceContainer;

	typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric, KernelFunction, DTL2> VCMOnSurface;
	typedef KSpace::Surfel Surfel;
	typedef VoronoiMap<Space, NotPointPredicate, Metric> VoronoiMap;
	typedef Eigen::MatrixXd MatrixXd;

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (corresponding volume)")
		("output,o", po::value<std::string>(), "output skeleton filename")
		("skeleton,s", po::value<std::string>(), "vol file (medial axis)")
		("delta,d", po::value<double>()->default_value(1), "delta for ball radius")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		("radiusInside,R", po::value<double>()->default_value(10), "radius of the ball inside voronoi cell")
		("radiusNeighbour,r", po::value<double>()->default_value(10), "radius of the ball for the neighbourhood")
		("angleThreshold,a", po::value<double>()->default_value(0.1), "anglem threshold")
		;

	bool parseOK=true;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, general_opt), vm);
	} catch(const std::exception& ex){
		parseOK=false;
		trace.info()<< "Error checking program options: "<< ex.what()<< endl;
	}
	po::notify(vm);
	if( !parseOK || vm.count("help")||argc<=1)
	{
		std::cout << "Usage: " << argv[0] << " [input]\n"
				  << "Display volume file as a voxel set by using QGLviewer"<< endl
				  << general_opt << "\n";
		return 0;
	}
	if(!vm.count("input"))
	{
		trace.error() << " The file name was not defined" << endl;
		return 0;
	}

	string outFilename = vm["output"].as<std::string>();
	string inputFilename = vm["input"].as<std::string>();
	string skeletonFilename = vm["skeleton"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	double R = vm["radiusInside"].as<double>();
	double r = vm["radiusNeighbour"].as<double>();
	double delta = vm["delta"].as<double>();
	double angleThreshold = vm["angleThreshold"].as<double>();
//	bool isDT = vm["skeleton"].as<bool>();
	bool isDT = true;

	QApplication application(argc,argv);
	Viewer3D<> viewer;
	viewer.show();

	Image volume = VolReader<Image>::importVol(inputFilename);
	Z3i::Domain domainVolume = volume.domain();
	Z3i::DigitalSet setVolume(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
												  thresholdMin-1, thresholdMax);
	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);


	Image skeleton = VolReader<Image>::importVol(skeletonFilename);
	Z3i::Domain domainSkeleton = skeleton.domain();
	Z3i::DigitalSet setSkeleton(domainSkeleton);
	Z3i::DigitalSet branchingPoints(domainSkeleton);
	SetFromImage<Z3i::DigitalSet>::append<Image>(setSkeleton, skeleton, thresholdMin-1, thresholdMax);

	Z3i::DigitalSet existingSkeleton = CurveAnalyzer::ensureConnexity(setSkeleton);
	typedef StandardDSS6Computer<vector<Point>::iterator,int,8> SegmentComputer;
	typedef GreedySegmentation<SegmentComputer> Segmentation;
	vector<Point> existingSkeletonOrdered;
	Z3i::Point p = (*existingSkeleton.begin());
	// We have to visit the direct neighbours in order to have a container with voxels
	// ordered sequentially by their connexity
	// Otherwise we have a point container with points which are not neighbours
	// and this impairs maximal segment recognition
	typedef MetricAdjacency Graph;
	typedef DepthFirstVisitor<Graph, set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;
	typedef GraphVisitorRange<Visitor> VisitorRange;
	Graph graph;
	Visitor visitor( graph, p );
	MyNode node;

    unsigned int previous = 0;
	while ( !visitor.finished() )
	{
  		node = visitor.current();
		if ( existingSkeleton.find(node.first) != existingSkeleton.end() ) { //is inside domain
		    if (node.second <= previous) {
				vector<Z3i::Point> neighbors;
				back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
				MetricAdjacency::writeNeighbors(inserter, node.first);
				for (const Z3i::Point& n : neighbors) {
					if (find(existingSkeletonOrdered.begin(), existingSkeletonOrdered.end(), n) != existingSkeletonOrdered.end()) {
						branchingPoints.insert(n);
						existingSkeletonOrdered.push_back(n);
					}
				}
			}
			previous = node.second;
			existingSkeletonOrdered.push_back(node.first);
			visitor.expand();
		}
		else
			visitor.ignore();
	}
	SegmentComputer algo;
	Segmentation s(existingSkeletonOrdered.begin(), existingSkeletonOrdered.end(), algo);
	s.setMode("MostCentered++");
	typename Segmentation::SegmentComputerIterator itseg = s.begin();
	typename Segmentation::SegmentComputerIterator end = s.end();
//	Z3i::DigitalSet branchingPoints = branchingPointDetection(s, existingSkeletonOrdered, domainVolume);
//	Z3i::DigitalSet branchingPoints = detectCriticalPoints(existingSkeleton);

	vector<Z3i::DigitalSet> edgeGraph = CurveAnalyzer::constructGraph(existingSkeletonOrdered, branchingPoints);
	vector<Z3i::Point> endPoints = CurveAnalyzer::findEndPoints(existingSkeleton);
	vector<GraphEdge*> hierarchicalGraph = CurveAnalyzer::hierarchicalDecomposition(edgeGraph, endPoints, branchingPoints);


	//Display points
//	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << branchingPoints;
	// for (const Z3i::Point& p : endPoints) {
	// 	viewer << CustomColors3D(Color::Green, Color::Green) << p;
	// }
//	 viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	// for (GraphEdge* edge : hierarchicalGraph) {
	// 	int label = edge->getLabel();
	// 	int r = (label * 64) % 256;
	// 	int g = (label* 128)%256;
	// 	int b = (label* 192)%256;
	// 	viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << edge->pointSet();
	// }
	// viewer << Viewer3D<>::updateDisplay;
	// application.exec();
	// return 0;


	set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount>> setVolumeWeighted;
	Image volumeBinary(volume.domain());
	for (auto it = volume.domain().begin(), ite = volume.domain().end(); it != ite; ++it) {
		if (volume(*it) >= thresholdMin && volume(*it) <= thresholdMax)
			volumeBinary.setValue(*it, 255);
		else
			volumeBinary.setValue(*it, 0);
	}

	vector<Point> vPoints;
	Z3i::DigitalSet skeletonPoints(domainVolume);
	ThresholdedImage binarizer(volume, thresholdMin-1, thresholdMax);
	BackgroundPredicate backgroundPredicate(binarizer);
	DTL2 dt(&volume.domain(), &binarizer, &Z3i::l2Metric);
	DTL2::Domain domainDT = dt.domain();
	for (auto it = domainDT.begin(), ite = domainDT.end(); it != ite; ++it) {
		double value = dt(*it);
		if (value > 0) {
			setVolumeWeighted.insert(new WeightedPointCount(*it, value));
			checkPointForMedialAxis(dt, vPoints, *it);
		}
	}

	vector<Z3i::DigitalSet> planes;
	for (const Z3i::DigitalSet& edge : edgeGraph) {
		if (edge.size() == 0) continue;
		vector<pair<Z3i::Point, double>> pointToAreas = areaProfile<VCM, KernelFunction>(edge, setVolumeWeighted);
		Z3i::Point point = pointsVaryingNeighborhood(pointToAreas);
		if (point == Z3i::Point()) continue;
		Z3i::DigitalSet plane = associatedPlane<VCM, KernelFunction>(point, edge, domainVolume, setVolumeWeighted);
		planes.push_back(plane);
	}
	vector<Z3i::DigitalSet> subVolumes = computeSubVolumes(setVolume, planes);
	trace.info() << subVolumes.size() << endl;
	for (auto it = subVolumes.begin(), ite = subVolumes.end(); it != ite; ++it) {
		int r = rand() % 256, g = rand() % 256, b = rand() % 256;
		viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << *it;
	}

	for (const Z3i::DigitalSet& sub : subVolumes) {
		vector<Z3i::DigitalSet> delineating = delineatingPlanes(sub, planes);
		vector<Z3i::DigitalSet> pairPlanes = twoClosestPlanesHausdorff(delineating);
		trace.info() << pairPlanes.size() << endl;
		for (const Z3i::DigitalSet& plane : pairPlanes) {
			viewer << CustomColors3D(Color::Black, Color::Black) << plane;
		}
	}
	// viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	// viewer << CustomColors3D(Color(0,0,120,20), Color(0,0,120,20)) << setVolume;
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;


	// for (const Z3i::DigitalSet& edge : edgeGraph) {
	// 	if (edge.size() == 0) continue;
	// 	map<Z3i::Point, double> pointToAreas = areaProfile<VCM, KernelFunction>(edge, setVolumeWeighted);
	// 	double minArea = min_element(pointToAreas.begin(), pointToAreas.end(), [&](const pair<Z3i::Point, double>& pair1, const pair<Z3i::Point, double>& pair2) {
	// 			return pair1.second < pair2.second;
	// 		})->second;
		
	// 	double maxArea = max_element(pointToAreas.begin(), pointToAreas.end(), [&](const pair<Z3i::Point, double>& pair1, const pair<Z3i::Point, double>& pair2) {
	// 			return pair1.second < pair2.second;
	// 		})->second;

	// 	GradientColorMap<double> gradient( minArea, maxArea);
	// 	gradient.addColor(DGtal::Color::Blue);
	// 	gradient.addColor(DGtal::Color::Green);
	// 	gradient.addColor(DGtal::Color::Yellow);
	// 	gradient.addColor(DGtal::Color::Red);
	// 	Z3i::Point point = pointsVaryingNeighborhood(pointToAreas);
	// 	pair<Z3i::Point, Z3i::RealPoint> pair2 = pointToNormal<VCM, KernelFunction>(point, edge, domainVolume, setVolumeWeighted);
	// 	Z3i::Point currentP = pair2.first;
	// 	Z3i::RealPoint normalP = pair2.second;
	// 	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << currentP;
	// 	vector<Z3i::RealPoint> plane = SliceUtils::computePlaneFromNormalVector(normalP, currentP);
	// 	viewer.setFillColor(Color::Red);
	// 	double factor = 15;
	// 	viewer.addQuad(currentP+(plane[0]-currentP)*factor, currentP+(plane[1]-currentP)*factor, currentP+(plane[2]-currentP)*factor, currentP+(plane[3]-currentP)*factor);
	// 	viewer << Viewer3D<>::updateDisplay;
	// 	qApp->processEvents();		 
	// }
	// viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	// viewer << CustomColors3D(Color(0,0,120,20), Color(0,0,120,20)) << setVolume;
	// viewer << Viewer3D<>::updateDisplay;
	// application.exec();
	// return 0;

	
	
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	SaddleComputer<DTL2, BackgroundPredicate> saddleComputer(setVolume, dt, backgroundPredicate, R, r, delta);
	Z3i::DigitalSet saddlePoints = saddleComputer.extractSaddlePoints(setVolume);
	vector<Z3i::Object26_6> objSaddle = saddleComputer.saddleConnectedComponents(saddlePoints);
	Z3i::DigitalSet maxCurvaturePoints = saddleComputer.saddlePointsToOnePoint<Matrix>(objSaddle);


	Z3i::Point dummy;

	//Z3i::DigitalSet saddleCurrentCC = s.pointSet();
	map<Z3i::Point, vector<pair<Z3i::Point, Z3i::RealPoint>>> saddleToPlane = constructSubVolumeWithPlanes<VCM, KernelFunction>(hierarchicalGraph, setVolumeWeighted, domainVolume, objSaddle);
	for (const auto& pairS : saddleToPlane) {
		vector<pair<Z3i::Point, Z3i::RealPoint>> mapPointToNormal = pairS.second;
		for (const auto & pair2 : mapPointToNormal) {
			Z3i::Point currentP = pair2.first;
			Z3i::RealPoint normalP = pair2.second;
			viewer << CustomColors3D(Color::Yellow, Color::Yellow) << currentP;
			vector<Z3i::RealPoint> plane = SliceUtils::computePlaneFromNormalVector(normalP, currentP);
			viewer.setFillColor(Color::Red);
			double factor = 15;
			viewer.addQuad(currentP+(plane[0]-currentP)*factor, currentP+(plane[1]-currentP)*factor, currentP+(plane[2]-currentP)*factor, currentP+(plane[3]-currentP)*factor);
			viewer << Viewer3D<>::updateDisplay;
			qApp->processEvents();
		}
//			viewer.addLine(currentP, currentP+(normalP*6));
	}

	viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	viewer << CustomColors3D(Color(0,0,120,20), Color(0,0,120,20)) << setVolume;
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;

	//Construct VCM surface
	Metric l2;
	VCM vcmSurface(R, ceil(r), l2, false);
	KernelFunction chiSurface(1.0, r);
	vcmSurface.init(setVolume.begin(), setVolume.end());

	KernelFunction chi( 1.0, r );

	int i = 0;

	Z3i::RealPoint normal(0,0,1);
	Z3i::RealPoint previousNormal=normal;

	Z3i::DigitalSet connectedComponent3D(domainVolume);
	Z3i::RealPoint realCenter;
	Z3i::Point centerOfMass;

	double distanceMax = dt(*max_element(existingSkeleton.begin(), existingSkeleton.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
				return dt(one) < dt(two);
			}));

	trace.beginBlock("Computing skeleton");

	for (GraphEdge* graphEdge : hierarchicalGraph) {
		Z3i::DigitalSet edge = graphEdge->pointSet();
		trace.progressBar(i, edgeGraph.size());
		if (edge.size() == 0) continue;
		VCM vcm(R, ceil(r), l2, false);
		vcm.init(edge.begin(), edge.end());
		bool add = true;
		//Main loop to compute skeleton (stop when no vol points left to process)
		int cpt = 0;
		Z3i::DigitalSet smoothedEdge(edge.domain());
		for (const Z3i::Point& s : edge)
		{

			double radius = r;

			//Distance transform value for VCM radius
			if (isDT) {
				radius = dt(s);
				radius += delta;
				if (radius > 0) {
					chi = KernelFunction( 1.0, radius);
				}
			}

			// Compute discrete plane
			radius = edge.size() * 0.4;
			set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount>> subVolumeWeighted;
			Z3i::DigitalSet subVolume = constructSubVolumeWithTangent<VCM, KernelFunction>(s,
																				   hierarchicalGraph,
																				   edge,
																				   setVolumeWeighted,
																				   branchingPoints,
																				   saddlePoints,
																				   saddleComputer.getSurface());

			if (subVolume.size() < setVolumeWeighted.size()) {
				for (auto it = subVolume.begin(), ite=  subVolume.end(); it != ite; ++it) {
					viewer << CustomColors3D(Color(0,120,0,120), Color(0,120,0,120)) << *it;
				}
			}
			connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, subVolumeWeighted,
																 s, normal,
																 0, radius, distanceMax, false);

			//Z3i::DigitalSet surfelSet = extractRelevantSurfacePointsForVolVCM(connectedComponent3D, dt, s);
			Z3i::DigitalSet surfelSet = extractSurfacePointsWithEigenVectors(connectedComponent3D, dt, s, dt(s));
		    vector<Z3i::Point> surfels(surfelSet.begin(), surfelSet.end());

			chiSurface = KernelFunction(1.0, radius);
			vcmSurface.setMySmallR(radius);


			double radiusSurface = dt(s) + delta;
			Z3i::RealPoint normalSurface = VCMUtil::computeNormalFromVCM(s, vcmSurface, chiSurface, 0, Z3i::RealVector(), surfels);
			double dot = std::abs(normalSurface.dot(normal));
			bool keep = (dot > angleThreshold);

			realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);
			centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
			//Center of mass computation
			if (realCenter != Z3i::RealPoint()) {
//				centerOfMass = findMaxDTInSet(connectedComponent3D, dt, s);

				//Pruning
				if (graphEdge->getLabel() == 1 && !keep) {
					//viewer << CustomColors3D(Color::Yellow, Color::Yellow) << s;
					viewer << CustomColors3D(Color::Red, Color::Red);
					viewer.addLine(s, s+normal*6);
					viewer << CustomColors3D(Color::Blue, Color::Blue);
					viewer.addLine(s, s+normalSurface*6);
					// viewer << surfelSet;
					cpt++;
				}
				else
					smoothedEdge.insert(centerOfMass);

			}
			else {
				viewer << CustomColors3D(Color::Red, Color::Red) << s;
			}
		}
		if (cpt >= 0.5 * edge.size())
			viewer << CustomColors3D(Color::Yellow, Color::Yellow) << edge;
		else
			viewer << CustomColors3D(Color::Red, Color::Red) << smoothedEdge;
		viewer << Viewer3D<>::updateDisplay;
		qApp->processEvents();

		i++;
	}

	trace.endBlock();
//	viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;

	//second pass
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
	    viewer << CustomColors3D(Color(0,0,120,30), Color(0,0,120,30)) << (*it)->myPoint;
	}

	Image outImage(volume.domain());
	DGtal::imageFromRangeAndValue(skeletonPoints.begin(), skeletonPoints.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
