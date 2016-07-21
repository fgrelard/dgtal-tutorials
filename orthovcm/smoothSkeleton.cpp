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
#include "clustering/kmeans.h"

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

Z3i::DigitalSet reduceClustersToCenters(const Z3i::DigitalSet& clusters, const Z3i::DigitalSet& skeleton) {
	Z3i::DigitalSet newClusters(clusters.domain());

	Z3i::Object26_6 obj(Z3i::dt26_6, clusters);
	vector<Z3i::Object26_6> objects;
	back_insert_iterator<vector<Z3i::Object26_6>> inserter(objects);
    obj.writeComponents(inserter);
	for (const Z3i::Object26_6& o : objects) {
		Z3i::DigitalSet currentPointSet = o.pointSet();
		if (currentPointSet.size() > 1) {
			Z3i::Point candidate;
			int previous=0;
			for (const Z3i::Point& b : currentPointSet) {
				Z3i::DigitalSet difference = skeleton;
				difference.erase(b);
				Z3i::Object26_6 differenceObj(Z3i::dt26_6, difference);
				vector<Z3i::Object26_6> objectsDifference;
				back_insert_iterator<vector<Z3i::Object26_6>> inserterDiff(objectsDifference);
				unsigned int nbCC = differenceObj.writeComponents(inserterDiff);
				if (nbCC > previous) {
					previous = nbCC;
					candidate = b;
				}
			}
			newClusters.insert(candidate);
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
																				 0, radius, radius*2, 6, false);
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

vector<pair<Z3i::Point, Z3i::RealPoint>> rotatedVectors(const pair<Z3i::Point, Z3i::RealPoint>& one,
												const pair<Z3i::Point, Z3i::RealPoint>& two) {
	vector<pair<Z3i::Point, Z3i::RealPoint>> rotated;
	Z3i::Point p1 = one.first;
	Z3i::RealPoint p2 =  two.first;
	Z3i::RealPoint n1 = one.second;
	Z3i::RealPoint n2 =  two.second;
	Z3i::RealPoint nrot = n1.crossProduct(n2);
	pair<Z3i::Point, Z3i::RealPoint> pair1 = make_pair(p1, n2.crossProduct(nrot));
	pair<Z3i::Point, Z3i::RealPoint> pair2 = make_pair(p2, n1.crossProduct(nrot));
	rotated.push_back(pair1);
	rotated.push_back(pair2);
	return rotated;
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
			// Z3i::RealPoint nrot2 = n2.crossProduct(n1);
			// double alpha = acos(n1.dot(n2));
			// float angleRotation = std::abs(M_PI/2 - alpha);
			// Eigen::Affine3f t(Eigen::AngleAxisf(angleRotation, Eigen::Vector3f(nrot[0], nrot[1], nrot[2])));
			// Eigen::Affine3f t2(Eigen::AngleAxisf(angleRotation, Eigen::Vector3f(nrot2[0], nrot2[1], nrot2[2])));
			// Eigen::Vector3f n1e(n1[0], n1[1], n1[2]);
			// Eigen::Vector3f n2e(n2[0], n2[1], n2[2]);
			// Eigen::Vector3f rotn1e = t.linear() * n1e;
			// Eigen::Vector3f rotn2e = t2.linear() * n2e;
			// Eigen::Vector3f otherRotn1e = t2.linear() * n1e;
			// Eigen::Vector3f otherRotn2e = t.linear() * n2e;
			// Z3i::RealPoint rotn1(rotn1e[0], rotn1e[1], rotn1e[2]);
			// Z3i::RealPoint rotn2(rotn2e[0], rotn2e[1], rotn2e[2]);
			// Z3i::RealPoint otherRotn1(otherRotn1e[0], otherRotn1e[1], otherRotn1e[2]);
			// Z3i::RealPoint otherRotn2(otherRotn2e[0], otherRotn2e[1], otherRotn2e[2]);
			// if (std::abs(otherRotn1.dot(otherRotn2)) < std::abs(rotn1.dot(rotn2))) {
			// 	rotn1 = otherRotn1;
			// 	rotn2 = otherRotn2;
			// }

			pointToVectors.clear();
			pointToVectors.push_back(make_pair(p1, n2.crossProduct(nrot)));
			pointToVectors.push_back(make_pair(p2, n1.crossProduct(nrot)));
			saddlesToPlanes[it->first] = pointToVectors;
		}
		else {
			saddlesToPlanes[it->first].clear();
		}
	}


	return saddlesToPlanes;
}

template <typename VCM, typename KernelFunction, typename Container, typename Domain>
vector<pair<Z3i::Point, double> > computeArea(const vector<Z3i::Point>& orientedEdge, const Container& setVolume, const Domain& domain) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;
	vector<pair<Z3i::Point, double>> areasEdge;
    auto itP = orientedEdge.begin();
	Z3i::Point start = *itP;
	Z3i::Point end = *(--orientedEdge.end());
	while (itP != orientedEdge.end()) {
		//if (*itP == start || *itP == end) {++itP;continue;} //End point effect
		Z3i::Point p = *itP;
		double radius = orientedEdge.size() * 0.4;
		VCM vcm(20, ceil(radius), l2, false);
		vcm.init(orientedEdge.begin(), orientedEdge.end());
		KernelFunction chi(1.0, radius);
		Z3i::RealPoint normal;
		Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domain, setVolume,
																			 p, normal,
																			 0, radius, radius*2, false);


		double area = connectedComponent3D.size();
		areasEdge.push_back(make_pair(p, area));
		++itP;
	}
	return areasEdge;
}


template <typename VCM, typename KernelFunction, typename Container>
vector<pair<Z3i::Point, double> > areaProfile(const Z3i::DigitalSet& setEdge, const Container& setVolume, const Z3i::Point& b) {
	auto domain = setEdge.domain();
	vector<Z3i::Point> orientedEdge = CurveAnalyzer::convertToOrientedEdge(setEdge, b);
	return computeArea<VCM, KernelFunction>(orientedEdge, setVolume, domain);
}




template <typename VCM, typename KernelFunction, typename Domain, typename WeightedContainer>
pair<Z3i::Point, Z3i::RealPoint> pointToNormal(const Z3i::Point&  point, const Z3i::DigitalSet& setVCM, const Domain& domain,
											  const WeightedContainer& setVolume) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space, 2> Metric;
	Metric l2;

	map<Z3i::Point, Z3i::RealPoint> aMap;
	double radius = setVCM.size()*0.4;
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

template <typename VCM, typename KernelFunction, typename Domain, typename Container, typename WeightedContainer>
Z3i::DigitalSet associatedPlane(const Z3i::Point& point, const Container& setVCM, const Domain& domain,
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

pair<Z3i::Point, double> pointsVaryingNeighborhood(const vector<pair<Z3i::Point, double>>& pointToAreas) {
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
				double factor = valueNeighbor / currentValue;
				if (factor > previousFactor // && factor > 2
					) {
					previousFactor = factor;
					candidate = p;
				}
			}
		}
	}
	return make_pair(candidate, previousFactor);
}

vector<pair<Z3i::Point, double> > areaVariationMeasure(const vector<pair<Z3i::Point, double>>& pointToAreas) {
	typedef MetricAdjacency<Z3i::Space, 3> MAdj;

	vector<pair<Z3i::Point, double>> areaVariation;
	for (const auto& pair : pointToAreas) {
		double previousFactor = 0;
		Z3i::Point p = pair.first;
		double currentValue = pair.second;
		double radius = sqrt(currentValue / M_PI);
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
	    MAdj::writeNeighbors(inserter, p);
		for (const Z3i::Point& n : neighbors) {
			auto nInMap = find_if(pointToAreas.begin(),pointToAreas.end(), [&](const std::pair<Z3i::Point, double>& pair) {
					return pair.first == n;
				});
			if (nInMap != pointToAreas.end()) {
				double valueNeighbor = nInMap->second;
				double factor = valueNeighbor / currentValue;
				if (factor > previousFactor // && factor > 2
					) {
					previousFactor = factor;
				}
			}
		}
		areaVariation.push_back(make_pair(p, previousFactor));
	}
	return areaVariation;
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
		Z3i::Point center = Statistics::extractCenterOfMass3D(plane);
		Ball<Z3i::Point> ball(center, 2);
		vector<Z3i::Point> neighbors = ball.pointsInBall();
		for (const Z3i::Point& n : neighbors) {
			dilatedPlane.insert(n);
		}
		// for (const Z3i::Point& pointPlane : plane) {
		// 	dilatedPlane.insert(pointPlane);
		// 	vector<Z3i::Point> neighbors;
		// 	back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
		// 	MetricAdjacency<Z3i::Space, 3>::writeNeighbors(inserter, pointPlane);
		// 	for (const Z3i::Point& n : neighbors) {
		// 		dilatedPlane.insert(n);
		// 	}
		// }
		for (const Z3i::Point& p : subVolume) {

			if (dilatedPlane.find(p) != dilatedPlane.end()) {
				delineating.push_back(plane);
				break;
			}
		}
	}
	return delineating;
}



vector<Z3i::DigitalSet> twoClosestPlanes(const vector<Z3i::DigitalSet>& delineatingPlanes) {
	double distanceMin = numeric_limits<double>::max();
	vector<Z3i::DigitalSet> pairPlanes;
	trace.info() << delineatingPlanes.size() << endl;
	for (size_t i = 0; i < delineatingPlanes.size(); i++) {
		for (size_t j = i+1; j < delineatingPlanes.size(); j++) {
			Z3i::DigitalSet first = delineatingPlanes[i];
			Z3i::DigitalSet second = delineatingPlanes[j];
			double dh = Distance::distanceSet(first, second);
			if (dh < distanceMin) {
				distanceMin = dh;
				vector<Z3i::DigitalSet> currentPair;
				currentPair.push_back(first);
				currentPair.push_back(second);
				pairPlanes = currentPair;
			}
		}
	}
	return pairPlanes;
}

pair<Z3i::Point, Z3i::Point> twoClosestPoints(const Z3i::DigitalSet& one, const Z3i::DigitalSet& two) {
	double distanceMax = std::numeric_limits<double>::max();
	Z3i::Point p1, p2;
	for (auto it = one.begin(), ite = one.end(); it != ite; ++it) {
		DGtal::Z3i::Point oneP = *it;
		DGtal::Z3i::Point closestPointInTheoretical =  *min_element(two.begin(), two.end(), [&](const DGtal::Z3i::Point& one, const DGtal::Z3i::Point& two) {
				return DGtal::Z3i::l2Metric(one, oneP) < DGtal::Z3i::l2Metric(two, oneP);
			});
		double distance = DGtal::Z3i::l2Metric(closestPointInTheoretical, oneP);
		if (distance < distanceMax) {
			distanceMax = distance;
			p1 = oneP;
			p2 = closestPointInTheoretical;
		}
	}
	return make_pair(p1, p2);

}

double gaussianDiff(double x, double sigma) {
	double diff = (-x/(pow(sigma, 3)*sqrt(2*M_PI)) * exp((-pow(x,2)/(2*pow(sigma, 2)))));
	return diff;
}



template <typename VCM, typename KernelFunction, typename Container, typename WeightedContainer, typename Domain>
void areaProfileToStdout(const vector<Container>& edgeGraph,
						 const WeightedContainer& setVolumeWeighted,
						 const Domain& domain) {
	int loop = 0;
	double sigma = 1;
	double w = 3*sigma;
	for (const Container& edge : edgeGraph) {
		if (edge.size() < 10) continue;
		vector<pair<Z3i::Point, double>> pointToAreas = computeArea<VCM, KernelFunction>(edge, setVolumeWeighted, domain);
		int j = 0;
		ofstream myfile;
		string outputname = "histos/histo" + to_string(loop) + ".csv";
		myfile.open (outputname);
		for (int i = 0+w; i < pointToAreas.size()-w; i++) {
			double area = pointToAreas[i].second;
			double sumDiffGaussian = 0;
			double mean = 0;
			for (int k = -w; k <= w; k++) {
				double currentArea = pointToAreas[i-k].second;
				mean += currentArea;
				sumDiffGaussian += currentArea * gaussianDiff(k, sigma);
			}
			mean /= (2 * w + 1);
			myfile << j << " " << area << " " << mean-area << " " << sumDiffGaussian << endl;
			j++;
		}
		myfile.close();
		loop++;
	}
}

vector<Z3i::DigitalSet> adjacentEdgesToBranchPoint(const vector<Z3i::DigitalSet>& edges,
												   const Z3i::Point& branchPoint,
												   const Z3i::DigitalSet& skeleton) {
	vector<Z3i::DigitalSet> adjacent;
	Z3i::Object26_6 objSkel(Z3i::dt26_6, skeleton);
	vector<Z3i::Point> neigh;
	back_insert_iterator<vector<Z3i::Point>> inserter(neigh);
	objSkel.writeNeighbors(inserter, branchPoint);

	for (const Z3i::DigitalSet& edge : edges) {
		for (const Z3i::Point& n : neigh) {
			if (edge.find(n) != edge.end()) {
				adjacent.push_back(edge);
				break;
			}
		}
	}
	return adjacent;
}


Z3i::DigitalSet restrictByDistanceToPoint(const Z3i::DigitalSet& original,
										  const Z3i::Point& p,
										  double distance) {
	Z3i::DigitalSet restricted(original.domain());
	for (const Z3i::Point& o : original) {
		if (Z3i::l2Metric(p, o) <= distance) {
			restricted.insert(o);
		}
	}
	return restricted;
}

pair<Z3i::Point, Z3i::RealPoint> intersectionBetweenPlanes(const pair<Z3i::Point, Z3i::RealPoint>& plane1,
														   const pair<Z3i::Point, Z3i::RealPoint>& plane2) {
	Z3i::RealPoint n1 = plane1.second;
	Z3i::RealPoint n2 = plane2.second;
	Z3i::Point p1 = plane1.first;
	Z3i::Point p2 = plane2.second;
	double d1 = n1[0] * p1[0] + n1[1] * p1[1] + n1[2] * p1[2];
	double d2 =	n2[0] * p2[0] + n2[1] * p2[1] + n2[2] * p2[2];
	Z3i::RealPoint direction = n1.crossProduct(n2);

	Eigen::Matrix2f A;
	Eigen::Vector2f b;

	A << n1[1], n1[2], n2[1], n2[2];
	b << d1, d2;

	Eigen::Vector2f x = A.colPivHouseholderQr().solve(b);
	Z3i::RealPoint solution(0, x[0], x[1]);
	return make_pair(solution, direction);
}

Z3i::DigitalSet createVolumeAroundPoint(const Z3i::DigitalSet& setVolume,
										const Z3i::Point& point, double radius) {
	Ball<Z3i::Point> ball(point, radius);
	vector<Z3i::Point> subVolumeInter = ball.intersection(setVolume);
	Z3i::DigitalSet subVolume(setVolume.domain());
	subVolume.insert(subVolumeInter.begin(), subVolumeInter.end());
	Z3i::Object26_6 obj(Z3i::dt26_6, subVolume);
	vector<Z3i::Object26_6> objects;
	back_insert_iterator<vector<Z3i::Object26_6>> inserter(objects);
	unsigned int nbCC = obj.writeComponents(inserter);
	if (nbCC > 1) {
		for (const Z3i::Object26_6& obj : objects) {
			Z3i::DigitalSet setObj = obj.pointSet();
			if (setObj.find(point) != setObj.end())
				return setObj;
		}
	}
	return subVolume;
}


Z3i::DigitalSet createSubVolume(const Z3i::DigitalSet& restrictedVolume,
								const Z3i::RealPoint& normal,
								const Z3i::Point& point) {
	Z3i::DigitalSet subVolume(restrictedVolume.domain());
	for (const Z3i::Point& p : restrictedVolume) {
		if (VCMUtil::abovePlane(p, normal, point))
			subVolume.insert(p);
	}
	return subVolume;
}

template <typename VCM, typename KernelFunction, typename Container>
Z3i::DigitalSet smoothedSkeletonPoints(const Z3i::DigitalSet& subVolume,
									   const Container& existingSkeleton,
									   const Z3i::DigitalSet& computedSkeleton) {

	typedef WeightedPointCount<Z3i::Point> WeightedPointCount;
	Z3i::DigitalSet smoothSkeleton(subVolume.domain());
	if (existingSkeleton.size() <= 2) {
		smoothSkeleton.insert(existingSkeleton.begin(), existingSkeleton.end());
		return smoothSkeleton;
	}

	set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount>> subVolumeWeighted;
	for (const Z3i::Point& subPoint : subVolume) {
		subVolumeWeighted.insert(new WeightedPointCount(subPoint, 1));
	}
	for (const Z3i::Point& cp : existingSkeleton) {
	    Z3i::DigitalSet connectedComponent3D = associatedPlane<VCM, KernelFunction>(cp, existingSkeleton, subVolume.domain(), subVolumeWeighted);
		Z3i::RealPoint realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);
		Z3i::Point centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
		bool stop = false;
		for (const Z3i::Point& p : computedSkeleton)
			if (Z3i::l2Metric(centerOfMass, p) <= sqrt(3))
				stop = true;
		if (stop) break;
		if (realCenter != Z3i::RealPoint()) {
			smoothSkeleton.insert(centerOfMass);
		}
	}
	return smoothSkeleton;
}

GraphEdge* firstEdge(const vector<GraphEdge*>& edgeGraph) {
	size_t maxSize = 0;
	GraphEdge* maxEdge;
	for (GraphEdge* edge : edgeGraph) {
		size_t currentSize = edge->pointSet().size();
		if (currentSize > maxSize) {
			maxSize = currentSize;
			maxEdge = edge;
		}
	}
	return maxEdge;
}

vector<vector<Z3i::Point> > sortedEdges(const vector<GraphEdge*>& unsortedEdges,
										GraphEdge* firstEdge,
										const Z3i::DigitalSet& branchingPoints) {

	vector<vector<Z3i::Point> > sortedEdges;
	Z3i::DigitalSet setEdge = firstEdge->pointSet();
	vector<Z3i::Point> endPoints = CurveAnalyzer::findEndPoints(setEdge);
	Z3i::Point start;
	for (const Z3i::Point& e : endPoints) {
		if (branchingPoints.find(e) == branchingPoints.end()) {
			start = e;
		}
	}
	queue<Z3i::Point> currentBranchingPoints;
	currentBranchingPoints.push(start);

	vector<GraphEdge*> neighbors;
	while (sortedEdges.size() < unsortedEdges.size()) {
		Z3i::Point b = currentBranchingPoints.front();
		currentBranchingPoints.pop();

		for (GraphEdge* edge : unsortedEdges) {
			if (find(neighbors.begin(), neighbors.end(), edge) != neighbors.end()) continue;
			Z3i::DigitalSet currentSet = edge->pointSet();
			if (currentSet.find(b) != currentSet.end()) {
				vector<Z3i::Point> orientedEdge = CurveAnalyzer::convertToOrientedEdge(currentSet, b);
				sortedEdges.push_back(orientedEdge);
				neighbors.push_back(edge);
				Z3i::Point end = *(--orientedEdge.end());
				for (GraphEdge* edge : unsortedEdges) {
					Z3i::DigitalSet currentSet = edge->pointSet();
					if (currentSet.find(end) != currentSet.end()) {
						currentBranchingPoints.push(end);
					}
				}
			}
		}
	}
	return sortedEdges;
}

bool planesIntersect(const Z3i::DigitalSet& plane,
					 const Z3i::DigitalSet& plane2) {
	for (const Z3i::Point& p : plane) {
		if (plane2.find(p) != plane2.end())
			return true;
	}
	return false;
}

template <typename VCM, typename KernelFunction, typename Domain, typename WeightedContainer>
std::vector<Z3i::DigitalSet> planesAlongEdge(const vector<Z3i::Point>& orientedEdge,
											 const Z3i::DigitalSet& edge,
											 const Domain& domain,
											 const WeightedContainer& setVolumeWeighted) {
	vector<Z3i::DigitalSet> planesEdge;
	for (const Z3i::Point& p : orientedEdge) {
		Z3i::DigitalSet plane = associatedPlane<VCM, KernelFunction> (p, edge, domain, setVolumeWeighted);
		planesEdge.push_back(plane);
	}
	return planesEdge;
}



void twoClosestPlanesNonIntersecting(Z3i::DigitalSet& plane,
									 Z3i::DigitalSet& plane2,
									 int& index,
									 int& index2,
									 const vector<Z3i::DigitalSet>& planesEdge,
									 const vector<Z3i::DigitalSet>& planesEdge2)
{
	for (int i = 0; i < planesEdge.size(); i ++) {
		Z3i::DigitalSet firstPlane = planesEdge[i];
		for (int j = 0; j < planesEdge2.size(); j++) {
			Z3i::DigitalSet secondPlane = planesEdge2[j];
			bool intersecting = planesIntersect(firstPlane, secondPlane);
			if (!intersecting ) {
				bool underIntersecting = false;

				for (int firstI = i; firstI < planesEdge.size(); firstI++) {
					Z3i::DigitalSet planeFirstUnder = planesEdge[firstI];
					for (int secondJ = j; secondJ < planesEdge2.size(); secondJ++) {
						Z3i::DigitalSet planeSecondUnder = planesEdge2[secondJ];
						underIntersecting |= planesIntersect(planeFirstUnder, planeSecondUnder);
					}
				}
				if (!underIntersecting) {
					index = i;
					index2 = j;
					plane = firstPlane;
					plane2 = secondPlane;
					return;
				}
			}
		}
	}
}

void twoClosestPlanesNonIntersecting(Z3i::DigitalSet& plane,
									 int& index,
									 const vector<Z3i::DigitalSet>& planesEdge,
									 const std::initializer_list<vector<Z3i::DigitalSet>>& otherPlanes)
{
	vector<size_t> doneIndexes;
	for (int i = 0; i < planesEdge.size(); i ++) {
		Z3i::DigitalSet firstPlane = planesEdge[i];
		for (int otherIndex = 0; otherIndex < otherPlanes.size(); otherIndex++) {
		    vector<Z3i::DigitalSet> otherPlanesEdge = otherPlanes.begin()[otherIndex];
			for (int j = 0; j < otherPlanesEdge.size(); j++) {
				Z3i::DigitalSet otherPlane = otherPlanesEdge[j];
				bool intersecting = planesIntersect(firstPlane, otherPlane);
				if (!intersecting) {
					bool underIntersecting = false;

					for (int firstI = i; firstI < planesEdge.size(); firstI++) {
						Z3i::DigitalSet planeFirstUnder = planesEdge[firstI];
						for (const vector<Z3i::DigitalSet>& restrictedPlanesEdge : otherPlanes) {
							for (int secondJ = j; secondJ < restrictedPlanesEdge.size(); secondJ++) {
								Z3i::DigitalSet planeSecondUnder = restrictedPlanesEdge[secondJ];
								underIntersecting |= planesIntersect(planeFirstUnder, planeSecondUnder);
							}
						}
					}
					if ((find(doneIndexes.begin(), doneIndexes.end(), otherIndex) == doneIndexes.end()) &&
						!underIntersecting) {
						doneIndexes.push_back(otherIndex);
						index = i;
						plane = firstPlane;
						if (doneIndexes.size() == otherPlanes.size()) return;
						break;
					}
				}
			}
		}

	}
}



///////////////////////////////////////////////////////////////////////////////
int main( int  argc, char**  argv )
{

	srand(time(NULL));
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
//	Z3i::DigitalSet branchingPoints(domainSkeleton);
	SetFromImage<Z3i::DigitalSet>::append<Image>(setSkeleton, skeleton, thresholdMin-1, thresholdMax);

	Z3i::DigitalSet existingSkeleton = CurveAnalyzer::ensureConnexity(setSkeleton);
	typedef StandardDSS6Computer<vector<Point>::iterator,int,8> SegmentComputer;
	typedef GreedySegmentation<SegmentComputer> Segmentation;
	Metric l2;
	//Z3i::DigitalSet branchingPoints = CurveAnalyzer::detectCriticalPoints(existingSkeleton);


	vector<Z3i::Point> endPointsV = CurveAnalyzer::findEndPoints(existingSkeleton);
	Z3i::Point p = (*endPointsV.begin());

	Z3i::DigitalSet branchingPoints(domainVolume);
	vector<Point> existingSkeletonOrdered = CurveAnalyzer::curveTraversalForGraphDecomposition(branchingPoints,
																							   existingSkeleton,
																							   p);


    //branchingPoints = CurveAnalyzer::detectCriticalPoints(existingSkeleton);
	vector<Z3i::Point> endPoints;
	for (const Z3i::Point& p : endPointsV) {
		if (branchingPoints.find(p) == branchingPoints.end())
			endPoints.push_back(p);
	}
	vector<Z3i::DigitalSet> edgeGraph = CurveAnalyzer::constructGraph(existingSkeletonOrdered, branchingPoints);
	vector<GraphEdge*> hierarchicalGraph = CurveAnalyzer::hierarchicalDecomposition(edgeGraph, endPoints, branchingPoints);



	//Display points
	// viewer << CustomColors3D(Color::Yellow, Color::Yellow) << branchingPoints;
	// for (const Z3i::Point& p : endPoints) {
	// 	viewer << CustomColors3D(Color::Green, Color::Green) << p;
	// }
	// viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	// for (GraphEdge* edge : hierarchicalGraph) {
	// 	int label = edge->getLabel();
	// 	int r = (label * 64) % 256;
	// 	int g = (label* 128)%256;
	// 	int b = (label* 192)%256;
	// 	viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << edge->pointSet();
	// }

	// for (const Z3i::DigitalSet& edge : edgeGraph) {
	// 	int r = rand() % 256, g = rand() % 256, b = rand() % 256;
	// 	viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << edge;
	// }
	// for (GraphEdge* edge : hierarchicalGraph) {
	// 	if (edge->getLabel() == 1) {
	//  		viewer << CustomColors3D(Color::Yellow, Color::Yellow);
	//  	}
	//  	else
	//  		viewer << CustomColors3D(Color::Red, Color::Red);
	//  	viewer << edge->pointSet();
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

	// GraphEdge* firstE = firstEdge(hierarchicalGraph);
	// vector<vector<Z3i::Point>> orderedEdges = sortedEdges (hierarchicalGraph, firstE, branchingPoints);
	// areaProfileToStdout<VCM, KernelFunction>(orderedEdges, setVolumeWeighted, domainVolume);
	// for (const vector<Z3i::Point>& edge : orderedEdges) {
	// 	int i = 0;
	// 	for (const Z3i::Point& p : edge) {
	// 		int r = i * 255. / edge.size();
	// 		viewer << CustomColors3D(Color(r,0,0,255), Color(r,0,0,255)) << p;
	// 		i++;
	// 	}
	// }
	// viewer << Viewer3D<>::updateDisplay;
	// application.exec();
	// return 0;


	vector<Z3i::DigitalSet> planes;
	/* Working with 3 planes */
	Z3i::DigitalSet processedEdges(existingSkeleton.domain());
	Z3i::DigitalSet skeletonPoints(setVolume.domain());
	for (const Z3i::Point& b : branchingPoints) {
		vector<Z3i::DigitalSet> adjacentEdges = adjacentEdgesToBranchPoint(edgeGraph, b, existingSkeleton);
		vector<Z3i::DigitalSet> junctions;
		vector<Z3i::DigitalSet> restrictedAdjacentEdges;
		double radius = max_element(adjacentEdges.begin(), adjacentEdges.end(), [&](const Z3i::DigitalSet& e1,
		 																			const Z3i::DigitalSet& e2) {
										return e1.size() < e2.size();
									})->size()*0.5;
		for (const Z3i::DigitalSet& adjacentE : adjacentEdges) {
			double radius = adjacentE.size()*0.8;
			DGtal::Z3i::DigitalSet restrictedAdj = restrictByDistanceToPoint(adjacentE, b, radius);
			restrictedAdjacentEdges.push_back(restrictedAdj);
		}
		if (adjacentEdges.size() > 3  || dt(b) <= sqrt(3)
			) {
			for (const Z3i::DigitalSet& restrictedAdj : restrictedAdjacentEdges) {
				processedEdges.insert(restrictedAdj.begin(), restrictedAdj.end());
			}
			for (const Z3i::DigitalSet& restrict : restrictedAdjacentEdges)
				viewer << CustomColors3D(Color::Red, Color::Red) << restrict;

			continue;
		}
		if (adjacentEdges.size() == 3) {
			vector<pair<Z3i::Point, double>> pointsVarying;
			vector<pair<Z3i::Point, Z3i::RealPoint>> pointsToNormals;
			vector<vector< pair< Z3i::Point, double > > > areasAdjacentEdges;
			double referenceArea = std::numeric_limits<double>::max();
			int index = -1;
			vector<Z3i::DigitalSet> planes;
			vector<vector<vector<Z3i::Point> > > clusters;
			vector<vector<vector<double> > > clusterValues;
			vector<Z3i::Point> cuttingPoints;
			for (int i = 0, end = restrictedAdjacentEdges.size(); i < end; i++) {
				Z3i::DigitalSet restrictEdge = restrictedAdjacentEdges[i];
				if (restrictEdge.size() == 0) continue;
				vector< pair< Z3i::Point, double > > pointToAreas = areaProfile< VCM, KernelFunction > (restrictEdge, setVolumeWeighted, b);
				if (pointToAreas.size() == 0) continue;
				areasAdjacentEdges.push_back(pointToAreas);
				pair<Z3i::Point, double> cuttingPoint = pointsVaryingNeighborhood (pointToAreas);
				cuttingPoints.push_back(cuttingPoint.first);
				vector<pair<Z3i::Point, double> > areaVariation = areaVariationMeasure(pointToAreas);

				vector<double> areaValues;
				for (const pair<Z3i::Point, double> pair : areaVariation) {
					areaValues.push_back(pair.second);
				}
			    vector<vector<double> > currentCluster = KMeans::kmeansAlgorithm(areaValues, 2);
				clusterValues.push_back(currentCluster);
				Z3i::DigitalSet plane = associatedPlane<VCM, KernelFunction> (cuttingPoint.first, restrictEdge, domainVolume, setVolumeWeighted);
				planes.push_back(plane);

				vector<vector<Z3i::Point>> associatedPoints(2);
				for (const pair<Z3i::Point, double>& pair : areaVariation) {
					if (find(currentCluster[0].begin(), currentCluster[0].end(), pair.second) != currentCluster[0].end())
						associatedPoints[0].push_back(pair.first);
					else
						associatedPoints[1].push_back(pair.first);
				}
				clusters.push_back(associatedPoints);
				// double maxAreaVar  = max_element(areaVariation.begin(), areaVariation.end(), [&](const pair<Z3i::Point, double>& one,
				// 																				 const pair<Z3i::Point, double>& two) {
				// 									 return one.second < two.second;
				// 								 })->second;
				// HueShadeColorMap<double, 1> hueShade(0.0, maxAreaVar);
				// Z3i::Point candidate;
				// for (const pair<Z3i::Point, double>& pair : areaVariation) {
				// 	if (pair.second > 0.9*maxAreaVar)
				// 		candidate = pair.first;
				// }
				// viewer << CustomColors3D(Color::Yellow, Color::Yellow) << candidate;
				double areaB = pointToAreas[0].second;
				if (areaB < referenceArea) {
					referenceArea = areaB;
					index = i;
				}
			}
			int indexDiscard = -1;
			if (clusterValues.size() == 3) {
				double minDiff = numeric_limits<double>::max();
				for (int i = 0, end = clusterValues.size(); i < end; i++) {
					vector<vector<double> > currentCluster = clusterValues[i];
					if (currentCluster[0].size() ==0 || currentCluster[1].size() == 0) continue;
					double firstValue = currentCluster[0][0];
					double secondValue = currentCluster[1][0];
					double diff = abs(firstValue - secondValue);
					if (diff < minDiff) {
						minDiff = diff;
						indexDiscard = i;
					}
				}

				trace.info() << indexDiscard << endl;
				cuttingPoints.erase(cuttingPoints.begin() + indexDiscard);
				planes.erase(planes.begin() + indexDiscard);

				Z3i::Point candFirst =  cuttingPoints[0], candSecond = cuttingPoints[1];
				Z3i::DigitalSet planeFirst = planes[0];
				Z3i::DigitalSet planeSecond = planes[1];
				for (int i = 0, end = planes.size(); i < end; i++) {
					Z3i::DigitalSet firstPlane = planes[i];
					for (int j = i+1; j < planes.size(); j++) {
						if (i == indexDiscard || j == indexDiscard) continue;
						Z3i::DigitalSet secondPlane = planes[j];
						if (planesIntersect(firstPlane, secondPlane)) {
							vector<vector<Z3i::Point>> firstCluster = clusters[i];
							vector<vector<Z3i::Point>> secondCluster = clusters[j];
							Z3i::DigitalSet firstEdge = restrictedAdjacentEdges[i];
							Z3i::DigitalSet secondEdge = restrictedAdjacentEdges[j];
							std::vector<Z3i::Point> firstEdgeOriented = CurveAnalyzer::convertToOrientedEdge(firstEdge, b);
							std::vector<Z3i::Point> secondEdgeOriented = CurveAnalyzer::convertToOrientedEdge(secondEdge, b);

							for (const Z3i::Point& p  : firstEdgeOriented) {
								if (find(firstCluster[1].begin(), firstCluster[1].end(), p) != firstCluster[1].end()) {
									candFirst = p;
									planeFirst = associatedPlane<VCM, KernelFunction> (candFirst, firstEdge, domainVolume, setVolumeWeighted);
								}
							}
							for (const Z3i::Point& p  : secondEdgeOriented) {
								if (find(secondCluster[1].begin(), secondCluster[1].end(), p) != secondCluster[1].end()) {
									candSecond = p;
									planeSecond = associatedPlane<VCM, KernelFunction> (candSecond, secondEdge, domainVolume, setVolumeWeighted);
								}

							}
						}
					}
				}


				viewer << CustomColors3D(Color::Yellow, Color::Yellow) << candFirst << candSecond;
				viewer << CustomColors3D(Color::Red, Color::Red) << planeFirst << planeSecond;
			}
		}
	}


	// 		vector<Z3i::DigitalSet> edgeAnalysis;
	// 		vector<vector<Z3i::DigitalSet> > planesRestrictedEdges;
	// 		for (int i = 0, end = areasAdjacentEdges.size(); i < end; i++) {
	// 			Z3i::DigitalSet restrictEdge = restrictedAdjacentEdges[i];
	// 			if (restrictEdge.size() == 0 // || i == index
	// 				) continue;
	// 			vector< pair< Z3i::Point, double > > pointToAreas = areasAdjacentEdges[i];
	// 			pair<Z3i::Point, double> point = pointsVaryingNeighborhood(pointToAreas);
	// 			pair<Z3i::Point, Z3i::RealPoint> ptoNPlane = pointToNormal<VCM, KernelFunction>(point.first, restrictEdge, domainVolume, setVolumeWeighted);
	// 			pointsToNormals.push_back(ptoNPlane);
	// 			edgeAnalysis.push_back(restrictEdge);
	// 			pointsVarying.push_back(point);
	// 			std::vector<Z3i::Point> restrictEdgeOriented = CurveAnalyzer::convertToOrientedEdge(restrictEdge, b);

	// 			std::vector<Z3i::DigitalSet> planesEdge = planesAlongEdge<VCM, KernelFunction> (restrictEdgeOriented, restrictEdge, domainVolume, setVolumeWeighted);
	// 			planesRestrictedEdges.push_back(planesEdge);
	// 		}


	// 		Z3i::DigitalSet cuttingPlane(domainVolume), cuttingPlane2(domainVolume), cuttingPlane3(domainVolume);
	// 		int firstIndex=0, firstIndex2 = 0, secondIndex = 0, secondIndex2 = 0, sumFirst = 0, sumSecond = 0, thirdIndex = 0, sumThird = 0;
	// 		int i  = 0;
	// 		int indexPlane1 = 0, indexPlane2 = 1;

	// 		vector<Z3i::DigitalSet> planesEdge = planesRestrictedEdges[0];
	// 		vector<Z3i::DigitalSet> planesEdge2 = planesRestrictedEdges[1];
	// 		vector<Z3i::DigitalSet> planesEdge3 = planesRestrictedEdges[2];

	// 		while ((i== 0 || firstIndex != 0 || secondIndex != 0 || thirdIndex != 0) &&
	// 			   sumFirst < planesEdge.size()-1 &&
	// 			   sumSecond < planesEdge2.size()-1 &&
	// 			   sumThird < planesEdge3.size()-1) {
	// 			std::vector<Z3i::DigitalSet> restrictedPl(planesEdge.begin() + sumFirst, planesEdge.end());
	// 			twoClosestPlanesNonIntersecting (cuttingPlane, firstIndex,
	// 											 restrictedPl,
	// 											 {planesEdge2, planesEdge3});
	// 			twoClosestPlanesNonIntersecting (cuttingPlane2, secondIndex,
	// 											 std::vector<Z3i::DigitalSet>(planesEdge2.begin() + sumSecond, planesEdge2.end()),
	// 											 {planesEdge, planesEdge3});
	// 			twoClosestPlanesNonIntersecting (cuttingPlane3, thirdIndex,
	// 											 std::vector<Z3i::DigitalSet>(planesEdge3.begin() + sumThird, planesEdge3.end()),
	// 											 {planesEdge, planesEdge2});
	// 			sumFirst += firstIndex;
	// 			sumSecond += secondIndex;
	// 			sumThird += thirdIndex;
	// 			i++;
	// 		}



	// 		//if (sumFirst == 0 || sumSecond == 0) continue;
	// 		viewer << CustomColors3D(Color::Red, Color::Red) << cuttingPlane << cuttingPlane2 << cuttingPlane3;
	// 		continue;

	// 		if (pointsVarying.size() < 2) continue;
	// 		//Find two max :  two planes
	// 		pair<Z3i::Point, double> maxiVarying = pointsVarying[0];
	// 		pair<Z3i::Point, double> maxiVarying2 = pointsVarying[1];

	// 		// double angleMax = 0.0;
	// 		// for (size_t i = 0; i < pointsToNormals.size(); i++) {
	// 		// 	Z3i::RealPoint ni = pointsToNormals[i].second;
	// 		// 	for (size_t j = i+1; j < pointsToNormals.size(); j++) {
	// 		// 		Z3i::RealPoint nj = pointsToNormals[j].second;
	// 		// 		double angle = ni.cosineSimilarity(nj);
	// 		// 		double otherAngle = ni.cosineSimilarity(-nj);
	// 		// 		angle = (angle < otherAngle) ? angle : otherAngle;
	// 		// 		if (angle > angleMax) {
	// 		// 			angleMax = angle;
	// 		// 			maxiVarying.first = pointsToNormals[i].first;
	// 		// 			maxiVarying2.first = pointsToNormals[j].first;
	// 		// 		}
	// 		// 	}
	// 		// }
	// 		viewer << CustomColors3D(Color::Yellow, Color::Yellow) << maxiVarying2.first << maxiVarying.first;
	// 		if (maxiVarying.first == Z3i::Point() || maxiVarying2.first == Z3i::Point()) continue;


	// 		Z3i::DigitalSet restrictEdge = edgeAnalysis[0];
	// 		Z3i::DigitalSet restrictEdge2 = edgeAnalysis[1];
	// 		Z3i::DigitalSet restrictEdgeB(domainVolume);
	// 	    for (const Z3i::DigitalSet& rEdge : restrictedAdjacentEdges) {
	// 			if (!(CurveAnalyzer::sameSet(rEdge, restrictEdge) ||
	// 				  CurveAnalyzer::sameSet(rEdge, restrictEdge2)))
	// 				restrictEdgeB = rEdge;
	// 		}

	// 		std::vector<Z3i::Point> eEdge = CurveAnalyzer::findEndPoints(restrictEdge);
	// 		std::vector<Z3i::Point> eEdge2 = CurveAnalyzer::findEndPoints(restrictEdge2);

	// 		Z3i::Point candEdge = b;
	// 		for (const Z3i::Point& p : eEdge) {
	// 			if (branchingPoints.find(p) == branchingPoints.end())
	// 				candEdge = p;
	// 		}


	// 		Z3i::Point candEdge2 = b;
	// 		for (const Z3i::Point& p : eEdge2) {
	// 			if (branchingPoints.find(p) == branchingPoints.end())
	// 				candEdge2 = p;
	// 		}
	// 		std::vector<Z3i::Point> restrictEdgeOriented = CurveAnalyzer::convertToOrientedEdge(restrictEdge, candEdge);
	// 		std::vector<Z3i::Point> restrictEdge2Oriented = CurveAnalyzer::convertToOrientedEdge(restrictEdge2, candEdge2);
	// 		std::vector<Z3i::Point> restrictEdgeBOriented = CurveAnalyzer::convertToOrientedEdge(restrictEdgeB, b);


	// 		Z3i::DigitalSet plane = associatedPlane<VCM, KernelFunction>(maxiVarying.first, restrictEdge, domainVolume, setVolumeWeighted);
	// 		Z3i::DigitalSet plane2 = associatedPlane<VCM, KernelFunction>(maxiVarying2.first, restrictEdge2, domainVolume, setVolumeWeighted);
	// 		Z3i::DigitalSet planeToProject = associatedPlane<VCM, KernelFunction>(b, restrictEdgeB, domainVolume, setVolumeWeighted);
	// 		if (plane.size() == 0 || plane2.size() == 0)
	// 			continue;


	// 		pair<Z3i::Point, Z3i::Point> closestPointsInter = twoClosestPoints(plane, plane2);

	// 		Z3i::Point current = closestPointsInter.first;
	// 		Z3i::Point current2 = closestPointsInter.second;


	// 		pair<Z3i::Point, Z3i::RealPoint> ptoNPlane = pointToNormal<VCM, KernelFunction>(maxiVarying.first, restrictEdge, domainVolume, setVolumeWeighted);
	// 		pair<Z3i::Point, Z3i::RealPoint>  ptoNPlane2 = pointToNormal<VCM, KernelFunction>(maxiVarying2.first, restrictEdge2, domainVolume, setVolumeWeighted);
	// 		pair<Z3i::Point, Z3i::RealPoint>  ptoNPlaneToProject = pointToNormal<VCM, KernelFunction>(b, restrictEdge2, domainVolume, setVolumeWeighted);
	// 		Z3i::RealPoint normalRot = ptoNPlane.second.crossProduct(ptoNPlane2.second);
	// 		Z3i::RealPoint normalPlane = normalRot.crossProduct(ptoNPlane2.second);
	// 		Z3i::RealPoint normalPlane2 = normalRot.crossProduct(ptoNPlane.second);

	// 		Z3i::DigitalSet newPlane(domainVolume), newPlane2(domainVolume);
	// 		double d = normalPlane[0] * current[0] + normalPlane[1] * current[1] + normalPlane[2] * current[2];
	// 		double omega = std::max(std::abs(normalPlane[0]), std::max(std::abs(normalPlane[1]), std::abs(normalPlane[2])));
	// 		VCMUtil::extractConnectedComponent3D(newPlane, domainVolume, setVolumeWeighted, normalPlane, current, d, omega);

	// 		d = normalPlane2[0] * current2[0] + normalPlane2[1] * current2[1] + normalPlane2[2] * current2[2];
	// 		omega = std::max(std::abs(normalPlane2[0]), std::max(std::abs(normalPlane2[1]), std::abs(normalPlane2[2])));
	// 		VCMUtil::extractConnectedComponent3D(newPlane2, domainVolume, setVolumeWeighted, normalPlane2, current2, d, omega);


	// 		Z3i::RealVector delineatePlane = (b - current).getNormalized();
	// 		Z3i::RealVector delineatePlane2 = (b -current2).getNormalized();
	// 		normalPlane = (delineatePlane.dot(normalPlane) < 0) ? -normalPlane : normalPlane;
	// 		normalPlane2 = (delineatePlane2.dot(normalPlane2) < 0) ? -normalPlane2 : normalPlane2;

	// 		Z3i::DigitalSet delineatedNewPlane(newPlane.domain());
	// 		for (const Z3i::Point& p : newPlane) {
	// 			if (VCMUtil::abovePlane(p, normalPlane2, current))
	// 				delineatedNewPlane.insert(p);
	// 		}


	// 		Z3i::DigitalSet delineatedNewPlane2(newPlane2.domain());
	// 		for (const Z3i::Point& p : newPlane2) {
	// 			if (VCMUtil::abovePlane(p, normalPlane, current2))
	// 				delineatedNewPlane2.insert(p);
	// 		}

	// 		planes.push_back(newPlane);
	// 		planes.push_back(newPlane2);


	// 		{
	// 		Z3i::DigitalSet restrictedVolumePlane = createVolumeAroundPoint(setVolume, b, radius*1.5);
	// 		Z3i::DigitalSet subVolume =  createSubVolume (restrictedVolumePlane, normalPlane, current);
	// 		Z3i::DigitalSet subVolume2 = createSubVolume (restrictedVolumePlane, normalPlane2, current);
	// 		Z3i::DigitalSet subVolumeUnder = createSubVolume (restrictedVolumePlane, -normalPlane2, current);
	// 		Z3i::DigitalSet subVolumeUnder2 = createSubVolume (restrictedVolumePlane, -normalPlane, current);
	// 		subVolume.insert(subVolumeUnder.begin(), subVolumeUnder.end());
	// 		subVolume2.insert(subVolumeUnder2.begin(), subVolumeUnder2.end());

	// 		std::vector<Z3i::Point> edgesVolume;
	// 	    edgesVolume.insert(edgesVolume.end(), restrictEdge2Oriented.begin(), restrictEdge2Oriented.end());
	// 		edgesVolume.insert(edgesVolume.end(), restrictEdgeBOriented.begin(), restrictEdgeBOriented.end());

	// 		std::vector<Z3i::Point> edgesVolume2;
	// 	    edgesVolume2.insert(edgesVolume2.end(), restrictEdgeOriented.begin(), restrictEdgeOriented.end());
	// 		edgesVolume2.insert(edgesVolume2.end(), restrictEdgeBOriented.begin(), restrictEdgeBOriented.end());

	// 		Z3i::DigitalSet smoothedSkeleton = smoothedSkeletonPoints<VCM, KernelFunction> (subVolume, edgesVolume, Z3i::DigitalSet(setVolume.domain()));
	// 		Z3i::DigitalSet smoothedSkeleton2 = smoothedSkeletonPoints<VCM, KernelFunction> (subVolume2, edgesVolume2, smoothedSkeleton);

	// 		viewer << CustomColors3D(Color::Blue, Color::Blue) << smoothedSkeleton;
	// 		viewer << CustomColors3D(Color::Blue, Color::Blue) << smoothedSkeleton2;
	// 		skeletonPoints.insert(smoothedSkeleton.begin(), smoothedSkeleton.end());
	// 		skeletonPoints.insert(smoothedSkeleton2.begin(), smoothedSkeleton2.end());
	// 		viewer << CustomColors3D(Color::Red, Color::Red) << delineatedNewPlane << delineatedNewPlane2;
	// 		for (const Z3i::DigitalSet& restrictedAdj : restrictedAdjacentEdges)
	// 			processedEdges.insert(restrictedAdj.begin(), restrictedAdj.end());
	// 		}


	// 		// vector<Z3i::RealPoint> planeToDisplay = SliceUtils::computePlaneFromNormalVector(normalPlane, closestPointsInter.first);
	// 		// vector<Z3i::RealPoint> planeToDisplay2 = SliceUtils::computePlaneFromNormalVector(normalPlane2, closestPointsInter.second);
	//  		// viewer.setFillColor(Color::Red);
	// 		// double factor = 15;

	//  		// viewer.addQuad(closestPointsInter.first+(planeToDisplay[0]-closestPointsInter.first)*factor, closestPointsInter.first+(planeToDisplay[1]-closestPointsInter.first)*factor, closestPointsInter.first+(planeToDisplay[2]-closestPointsInter.first)*factor, closestPointsInter.first+(planeToDisplay[3]-closestPointsInter.first)*factor);

	// 		// viewer.addQuad(closestPointsInter.second+(planeToDisplay2[0]-closestPointsInter.second)*factor, closestPointsInter.second+(planeToDisplay2[1]-closestPointsInter.second)*factor, closestPointsInter.second+(planeToDisplay2[2]-closestPointsInter.second)*factor, closestPointsInter.second+(planeToDisplay2[3]-closestPointsInter.second)*factor);
	// 	}
	// }
	 viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	 viewer << CustomColors3D(Color(0,0,120,20), Color(0,0,120,20)) << setVolume;
	 viewer << Viewer3D<>::updateDisplay;
	 application.exec();
	 return 0;
	vector<Z3i::DigitalSet> parts = computeSubVolumes(setVolume, planes);
	for (const Z3i::DigitalSet& part : parts) {
		int r = rand()%256, g = rand() % 256, b = rand() % 256;
		viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << part;
	}
	Z3i::DigitalSet notProcessed(existingSkeleton.domain());
	for (const Z3i::Point& s : existingSkeleton) {
		if (processedEdges.find(s) == processedEdges.end())
			notProcessed.insert(s);
	}

	Z3i::Object26_6 objNotProcessed(Z3i::dt26_6, notProcessed);
	vector<Z3i::Object26_6> objects;
	back_insert_iterator<vector<Z3i::Object26_6> > inserter(objects);
	objNotProcessed.writeComponents(inserter);
	for (const Z3i::Object26_6& o : objects) {
		Z3i::DigitalSet currentSet = o.pointSet();
		Z3i::DigitalSet smoothedSkeleton = smoothedSkeletonPoints<VCM, KernelFunction>(setVolume, currentSet, Z3i::DigitalSet(setVolume.domain()));
		viewer << CustomColors3D(Color::Red, Color::Red) << smoothedSkeleton;
		skeletonPoints.insert(smoothedSkeleton.begin(), smoothedSkeleton.end());
	}
	Image outImage(volume.domain());
	DGtal::imageFromRangeAndValue(skeletonPoints.begin(), skeletonPoints.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);
	viewer << CustomColors3D(Color(0,0,120,20), Color(0,0,120,20)) << setVolume;
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;

}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
