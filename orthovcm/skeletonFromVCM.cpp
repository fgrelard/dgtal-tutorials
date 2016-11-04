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
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/geometry/curves/Naive3DDSSComputer.h"
#include "DGtal/geometry/curves/estimation/LambdaMST3D.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"

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
#include "geometry/SaddleComputer.h"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "geometry/DigitalPlane.h"
#include "geometry/DigitalPlaneSet.h"
#include "surface/SurfaceTraversal.h"
#include "geometry/CurveAnalyzer.h"
#include "Statistics.h"
#include "shapes/Ball.h"
///////////////////////////////////////////////////////////////////////////////
#include "ReverseHatPointFunction.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

class LabelledPoint : public WeightedPoint<Z3i::Point> {
	typedef WeightedPoint<Z3i::Point> Base;
	using Base::Base;
	friend bool operator<(const LabelledPoint& it, const LabelledPoint& other) {
		return (it.myPoint < other.myPoint);
	}
};


Z3i::DigitalSet computeShell(const Z3i::Point& center, const Z3i::DigitalSet& setVolume, double radiusInnerBall, double radiusOuterBall) {
	Ball<Z3i::Point> ballInnerBall(center, radiusInnerBall);
	std::vector<Z3i::Point> points = ballInnerBall.surfaceIntersection(setVolume);
	Z3i::DigitalSet shell(setVolume.domain());
	for (const Z3i::Point& p : points) {
		shell.insert(p);
	}
	return shell;
	// Z3i::DigitalSet shell(setVolume.domain());
	// for (auto it = setVolume.begin(), ite = setVolume.end();
	// 	 it != ite; ++it) {
	// 	Z3i::Point current = *it;
	// 	if (!(ballInnerBall.contains(current)) && ballOuterBall.contains(current)) {
	// 		shell.insert(current);
	// 	}
	// }

	// return shell;
}

Z3i::DigitalSet computeHalfShell(const Z3i::Point& center, const Z3i::RealVector& dirVector,
								 const Z3i::DigitalSet& plane,
								 const Z3i::DigitalSet& setVolume, double radiusInnerBall) {
	typedef BreadthFirstVisitor<Z3i::Object26_6, std::set<Z3i::Point> > Visitor;
	typedef Visitor::Node Node;

	Z3i::Object26_6 obj(Z3i::dt26_6, setVolume);
    Visitor visitor(obj, plane.begin(), plane.end());
	Z3i::DigitalSet shell(setVolume.domain());
	int radius = (int) radiusInnerBall;
	while (!visitor.finished()) {
		Node node = visitor.current();
		if (node.second > radius) break;
		if (node.second == radius && VCMUtil::abovePlane(node.first, dirVector, center))
			shell.insert(node.first);
		visitor.expand();

	}
    return shell;
}


unsigned int computeDegree(const Z3i::DigitalSet& shell) {
	typedef Z3i::Object26_6 ObjectType;

	ObjectType objectImage(Z3i::dt26_6, shell);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);
	unsigned int cpt = 0;

	for (const auto& obj : objects) {
		if (obj.size() > 1)
			cpt++;
	}
	return cpt;
}



template <typename Image>
double computeRadiusFromIntersection(const Image& volume, const Z3i::Point& point, const Z3i::RealPoint& normal,
									 double radius) {
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
//	Z2i::DigitalSet aSet = VCMUtil::extractConnectedComponent(processImage, center, 1, 255);
	Z2i::DigitalSet aSet(domainImage2D);
	for (const Z2i::Point& p : domainImage2D) {
		if (processImage(p) >= 1) {
			aSet.insert(p);
		}
	}
	Eigen::MatrixXd covmatrix = Statistics::computeCovarianceMatrix<Eigen::MatrixXd>(aSet);
	if (covmatrix.size() == 0) return 0;
	Z2i::RealVector projection = Statistics::extractEigenVector<Z2i::RealVector>(covmatrix, 1);
	Z2i::Point trackedPoint = PointUtil::trackPoint(center, aSet, projection);
	Z2i::Point otherTrackedPoint = PointUtil::trackPoint(center, aSet, -projection);
	double distance = Z2i::l2Metric(otherTrackedPoint, trackedPoint) / 2.0;
	return distance;
}


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



vector<Z3i::Point> findEndPoints(const Z3i::DigitalSet& set) {
	Z3i::Object26_6 objectSet(Z3i::dt26_6, set);
	vector<Z3i::Point> endPoints;
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
		objectSet.writeNeighbors(inserter, *it);
		if (neighbors.size() == 1)
			endPoints.push_back(*it);
	}
	return endPoints;
}



template <typename VCM, typename KernelFunction, typename Domain, typename Container>
pair<Z3i::DigitalSet, double> eigenValuesWithVCM(VCM& vcm, KernelFunction& chi, const Domain& domain,
												 const Z3i::Point& p, double radius,
												 const Container& setVolumeWeighted, const Z3i::DigitalSet& setVolume,
												 Viewer3D<>& viewer) {
	map<Z3i::DigitalSet, double> eigenValues;
	Z3i::RealPoint normal;
	double size = 0;
	Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domain, setVolumeWeighted,
																		 p, normal,
																		 0, radius, radius*10, 26, true);

	Ball<Z3i::Point> ball(p, radius);
	vector<Z3i::Point> vPoints = ball.intersection(setVolume);
	Z3i::RealVector eigenValue = VCMUtil::computeEigenValuesFromVCM(p, vcm, chi);
	double factor = (pow(vcm.R(), 4) * radius) / 8;
	double l0 = eigenValue[0];
	double l1 = eigenValue[1];
	double l2 = eigenValue[2] ;
	double eigenVar = l0 / (l0 + l1 + l2);
	double angle = M_PI/2 - vcm.vectorVariability(normal, chi, p);
	Z3i::DigitalSet discretePlane(connectedComponent3D.domain());
	for (const Z3i::Point& c : connectedComponent3D) {
		if (Z3i::l2Metric(p, c) <= radius)
			discretePlane.insert(c);
	}

	return make_pair(discretePlane, eigenVar);
}

template <typename DTL2>
void fillHoles(Z3i::DigitalSet& skeletonPoints, const DTL2& dt) {
	Z3i::Object26_6 objectImage(Z3i::dt26_6, skeletonPoints);
	vector<Z3i::Object26_6> skeletonCC;
	back_insert_iterator< std::vector<Z3i::Object26_6> > inserterCC( skeletonCC );
	objectImage.writeComponents(inserterCC);
	bool objectsToLink = true;
	while (objectsToLink) {
		objectsToLink = false;
		double distance = numeric_limits<double>::max();
		Z3i::Point currentLink, otherLink;
		int indexDelete1, indexDelete2;
		for (size_t i = 0; i < skeletonCC.size(); i++) {
			Z3i::DigitalSet currentCC = skeletonCC[i].pointSet();
			std::vector<Z3i::Point> endPointCurrent = CurveAnalyzer::findEndPoints(currentCC);
			for (size_t j = 0; j < skeletonCC.size(); j++) {
				if (i == j) continue;
				Z3i::DigitalSet otherCC = skeletonCC[j].pointSet();
				std::vector<Z3i::Point> endPointOther = CurveAnalyzer::findEndPoints(otherCC);
				for (auto itCurrent = endPointCurrent.begin(), itCurrentE = endPointCurrent.end();
					 itCurrent != itCurrentE; ++itCurrent) {
					Z3i::Point current = *itCurrent;
					for (auto itOther = endPointOther.begin(), itOtherE = endPointOther.end(); itOther != itOtherE;
						 ++itOther) {
						Z3i::Point other = *itOther;
						double currentDistance = Z3i::l2Metric(current, other);
						if (currentDistance > sqrt(3) &&
							currentDistance <= 2 * sqrt(3) &&
							currentDistance < distance) {
							distance = currentDistance;
							currentLink = current;
							otherLink = other;
							indexDelete1 = i;
							indexDelete2 = j;
							objectsToLink = true;
						}
					}
				}
			}
		}
		if (objectsToLink) {
			vector<Z3i::Point> link = PointUtil::linkTwoPoints26(currentLink, otherLink);
			Z3i::DigitalSet toKeep = skeletonCC[indexDelete1].pointSet();
			Z3i::DigitalSet toDelete = skeletonCC[indexDelete2].pointSet();
			toKeep.insert(link.begin(), link.end());
			toKeep.insert(toDelete.begin(),
						  toDelete.end());
			Z3i::Object26_6 toAdd(Z3i::dt26_6, toKeep);
			if (indexDelete1 == skeletonCC.size() - 1) {
				skeletonCC.pop_back();
				skeletonCC[indexDelete2] = skeletonCC[indexDelete1-1];
				skeletonCC.pop_back();
			}
			else if (indexDelete1 == skeletonCC.size() - 2) {
				skeletonCC[indexDelete2] = skeletonCC[indexDelete1+1];
				skeletonCC.pop_back();
				skeletonCC.pop_back();
			}
			else if (indexDelete2 == skeletonCC.size() - 1) {
				skeletonCC.pop_back();
				skeletonCC[indexDelete1] = skeletonCC[indexDelete2-1];
				skeletonCC.pop_back();
			}
		    else if (indexDelete2 == skeletonCC.size() - 2) {
				skeletonCC[indexDelete1] = skeletonCC[indexDelete2+1];
				skeletonCC.pop_back();
				skeletonCC.pop_back();
			}
			else {
				skeletonCC[indexDelete1] = skeletonCC[skeletonCC.size() - 1];
				skeletonCC[indexDelete2] = skeletonCC[skeletonCC.size() - 2];
				skeletonCC.pop_back();
				skeletonCC.pop_back();
			}
			skeletonCC.push_back(toAdd);
			bool isAdd = true;
			for (const auto & p: link) {
				if (dt(p) == 0)
					isAdd=false;
			}
			if (isAdd)
				skeletonPoints.insert(link.begin(), link.end());
		}
	}
}



Z3i::DigitalSet endPointCurves(const Z3i::DigitalSet& curves, const Z3i::DigitalSet& endPoints) {
	typedef Z3i::Object26_6 ObjectType;

	Z3i::DigitalSet endPointCurves(curves.domain());
	ObjectType objCurves(Z3i::dt26_6, curves);
	vector<ObjectType> curvesCC;
	back_insert_iterator<vector<ObjectType> > inserter(curvesCC);
	objCurves.writeComponents(inserter);
	for (const ObjectType& curve : curvesCC) {
		Z3i::DigitalSet curveSet = curve.pointSet();
		vector<Z3i::Point> localEndPoints = CurveAnalyzer::findEndPoints(curveSet);
		for (const Z3i::Point& e : localEndPoints) {
			if (endPoints.find(e) == endPoints.end() && curves.find(e) != curves.end())
				endPointCurves.insert(e);
		}
	}
	return endPointCurves;
}



Z3i::DigitalSet computeTraversedPoints(const Z3i::DigitalSet& setVolume,
									   const Z3i::Point& point,
									   const Z3i::RealPoint& normal ) {
	int scalar = 1;
	Z3i::DigitalSet traversed(setVolume.domain());

	Z3i::Point projection = point;
	while (setVolume.find(projection) != setVolume.end()) {
		projection = point + normal * scalar;
		traversed.insert(projection);
		scalar++;
	}
	return traversed;
}

template <typename DTL2, typename Container>
vector<Z3i::Point> linkWithBezierCurves(const map<Z3i::Point, Z3i::RealVector>& mapPointToNormal, Z3i::DigitalSet& skeletonPoints,
										const Container& setVolumeWeighted, const Z3i::DigitalSet& setVolume,
										const Z3i::Point& referencePoint,
										const Z3i::Point& otherPoint, const DTL2& dt) {
	Z3i::Domain domainVolume = skeletonPoints.domain();
	double radiusCurrentObject = dt(referencePoint) ;
	Z3i::RealPoint normalCurrentObject = mapPointToNormal.at(referencePoint);

	double radiusReference = dt(otherPoint);
	Z3i::RealPoint normalReference = (referencePoint - otherPoint).getNormalized();
	Z3i::RealPoint dirVectorCurrent = (otherPoint - referencePoint).getNormalized();
	if (normalCurrentObject.dot(dirVectorCurrent) < 0)
		normalCurrentObject = -normalCurrentObject;
	Z3i::RealPoint dirVectorReference = (referencePoint - otherPoint).getNormalized();
	if (normalReference.dot(dirVectorReference) < 0)
		normalReference = -normalReference;

	Z3i::DigitalSet traversedCurrent = computeTraversedPoints(setVolume, referencePoint, normalCurrentObject);
	Z3i::DigitalSet traversedReference = computeTraversedPoints(setVolume, otherPoint, normalReference);
	double distanceCR = numeric_limits<double>::max();
	Z3i::Point controlCurrent, controlReference;
	for (auto it = traversedCurrent.begin(), ite = traversedCurrent.end(); it != ite; ++it) {
		Z3i::Point nearest = extractNearestNeighborInSetFromPoint(traversedReference, *it);
		double currentDistance = Z3i::l2Metric(nearest, *it);
		if (currentDistance < distanceCR && currentDistance > sqrt(3)) {
			distanceCR = currentDistance;
			controlCurrent = *it;
			controlReference = nearest;
		}
	}
	vector<Z3i::Point> newPoints = PointUtil::bezierCurveDeCasteljau(otherPoint, referencePoint, controlReference, controlCurrent);
	return newPoints;
}




Z3i::DigitalSet shellPointsToShellAreas(const Z3i::DigitalSet& setVolume,
										const map<Z3i::Point, int>& shellPoints,
										const Z3i::DigitalSet& skeletonPoints) {
	Z3i::DigitalSet shellAreas(setVolume.domain());
	for (const Z3i::Point& p : setVolume) {
		Z3i::Point closestPoint = *min_element(skeletonPoints.begin(), skeletonPoints.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
					return Z3i::l2Metric(one, p) < Z3i::l2Metric(two, p);
			});
		if (shellPoints.find(closestPoint) != shellPoints.end())
			shellAreas.insert(p);
	}
	return shellAreas;
}




Z3i::DigitalSet pointsToLink(const Z3i::DigitalSet& setVolume,
							 const Z3i::DigitalSet& junction,
							 const Z3i::DigitalSet& endPoints,
							 int numberLimit) {
	typedef BreadthFirstVisitor<Z3i::Object26_6, std::set<Z3i::Point> > Visitor;
	typedef Visitor::Node Node;

	Z3i::DigitalSet pointsToLink(setVolume.domain());
	Z3i::Object26_6 obj(Z3i::dt26_6, setVolume);
    Visitor visitor(obj, junction.begin(), junction.end());
	Z3i::DigitalSet shell(setVolume.domain());
	int numberEndPointsFound = 0;
	while (!visitor.finished()) {
		if (numberEndPointsFound == numberLimit) break;
		Node node = visitor.current();
		if (endPoints.find(node.first) != endPoints.end()) {
			numberEndPointsFound++;
			pointsToLink.insert(node.first);
		}
		visitor.expand();
	}
	return pointsToLink;
}

template <typename DTL2>
Z3i::Point referencePointToLink(const Z3i::DigitalSet& setVolume,
								const Z3i::DigitalSet& pointsToLink,
								const DTL2& dt) {
	if (pointsToLink.size() <= 2 && pointsToLink.size() > 0) return (*pointsToLink.begin());
	double maxDistance = numeric_limits<double>::min();
	Z3i::Point reference;

	for (const Z3i::Point& p1 : pointsToLink) {
		double number = 0, sum = 0;
		for (const Z3i::Point& p2 : pointsToLink) {
			if (p1 == p2) continue;
			vector<Z3i::Point> points = PointUtil::linkTwoPoints(p1, p2);
			for (const Z3i::Point& l : points) {
				sum += dt(l);
			}
			number += points.size();
		}
		double sumDistanceDT = sum / number;
		if (sumDistanceDT > maxDistance) {
			maxDistance = sumDistanceDT;
			reference = p1;
		}
	}

	return reference;
}

Z3i::DigitalSet dilateJunctions(const Z3i::DigitalSet& setVolume,
								const Z3i::DigitalSet& skeletonPoints,
								const Z3i::DigitalSet& junctionAreas) {
	Z3i::DigitalSet dilatedJunction(junctionAreas.domain());
	Z3i::Object26_6 obj(Z3i::dt26_6, setVolume);
	for (const Z3i::Point& s : skeletonPoints) {
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point> > inserter(neighbors);
		obj.writeNeighbors(inserter, s);
		int cpt = 0;
		for (const Z3i::Point& n : neighbors) {
			if (junctionAreas.find(n) != junctionAreas.end())
				cpt++;
		}
		if (cpt > 0)
			dilatedJunction.insert(s);
	}
	return dilatedJunction;
}

map<Z3i::Point, int> shellPointsInArea(const map<Z3i::Point, int>& shellPoints,
									   const Z3i::DigitalSet& area) {
	map<Z3i::Point, int> shellPointsInArea;
	for (const pair<Z3i::Point, int>& sp : shellPoints) {
		if (area.find(sp.first) != area.end()) {
			shellPointsInArea.insert(sp);
		}
	}
	return shellPointsInArea;
}




Z3i::DigitalSet extractPart(const Z3i::Point& e,
									const Z3i::DigitalSet& endPoints,
									const vector<Z3i::DigitalSet>& parts) {
	Z3i::DigitalSet involvedPart(endPoints.domain());

	for (const Z3i::DigitalSet& part : parts) {
		if (part.find(e) != part.end())
			involvedPart = part;
	}
	return involvedPart;
}

Z3i::Point orientationNormal(const Z3i::DigitalSet& points,
							 const map<Z3i::Point, Z3i::RealVector> pointToNormal) {
	int referenceNumber = numeric_limits<int>::min();
	Z3i::Point cand = *(points.begin());
	for (const Z3i::Point& p : points) {
		Z3i::RealVector n = pointToNormal.at(p);
		int cpt = 0;
		for (const Z3i::Point& p2 : points) {
			if (p == p2) continue;
			Z3i::RealVector n2 = pointToNormal.at(p2);
			if (n.dot(n2) < 0)
				cpt++;
		}
		if (cpt > referenceNumber) {
			referenceNumber = cpt;
			cand = p;
		}
	}
	return cand;
}

void orientNormalsEndPoints(std::map<Z3i::Point, Z3i::RealVector>& mapPointToNormal,
							const Z3i::DigitalSet& endPoints,
							const vector<Z3i::DigitalSet>& parts) {
	for (const Z3i::Point& currentPoint : endPoints) {
		Z3i::DigitalSet involvedPart = extractPart(currentPoint, endPoints, parts);
		vector<Z3i::Point> endPointsParts = CurveAnalyzer::findEndPoints(involvedPart);

		auto iterator = mapPointToNormal.find(currentPoint);
		if (iterator == mapPointToNormal.end()) {
			continue;
		}
		Z3i::RealVector normal = mapPointToNormal.at(currentPoint);

		Z3i::Point otherPoint;
		for (const Z3i::Point& otherE : endPointsParts) {
			if (otherE != currentPoint)
				otherPoint = otherE;
		}
		if (otherPoint == Z3i::Point::zero) continue;
		Z3i::RealVector dirVector = (currentPoint - otherPoint).getNormalized();
		normal = (normal.dot(dirVector) < 0) ? -normal : normal;

		mapPointToNormal[currentPoint] = normal;
	}
}

vector<WeightedPointCount<Z3i::Point>* > distanceMapFromBranches(const Z3i::DigitalSet& skeletonPoints,
																 const vector<Z3i::DigitalSet>& branches,
																 const Z3i::DigitalSet& endPoints,
																 const vector<Z3i::DigitalSet>& parts,
																 const map<Z3i::Point, Z3i::RealPoint>& mapPointToNormal,
																 Viewer3D<>& viewer) {
	typedef WeightedPointCount<Z3i::Point> WPoint;
	typedef Naive3DDSSComputer < std::vector<Z3i::Point>::const_iterator, int, 8 > SegmentComputer;
	typedef SaturatedSegmentation<SegmentComputer> Segmentation;

	int label = 0;
	vector< WPoint* > distanceMap;
	for (const Z3i::Point& e : endPoints) {
		distanceMap.push_back(new WPoint(e, numeric_limits<double>::max(), label));
	}
	label++;

	Z3i::Object26_6 obj(Z3i::dt26_6, skeletonPoints);
	for (const Z3i::DigitalSet& junction : branches) {
		for (WPoint* dtp : distanceMap) {
			Z3i::Point currentPoint = dtp->myPoint;
			// Z3i::DigitalSet involvedPart = extractPart(currentPoint, endPoints, parts);
			// vector<Z3i::Point> endPointsParts = CurveAnalyzer::findEndPoints(involvedPart);

			auto iterator = mapPointToNormal.find(currentPoint);
			if (iterator == mapPointToNormal.end()) {
			    continue;
			}

			Z3i::RealVector normal = mapPointToNormal.at(currentPoint);

			double currentDistance = dtp->myWeight;
			DigitalPlane<Z3i::Space> digPlane(currentPoint, normal);
			Z3i::Point closest =  *min_element(junction.begin(), junction.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
					return Z3i::l2Metric(one, currentPoint) < Z3i::l2Metric(two, currentPoint);
			});
			bool add = false;
			for (const Z3i::Point& p : junction) {
				add |= digPlane.isPointAbove(p);
			}
			viewer << CustomColors3D(Color::Red, Color::Red);
			viewer.addLine(currentPoint, currentPoint+normal*10);
			double distance = Z3i::l2Metric(currentPoint, closest);
			if (currentDistance > distance && add) {
				dtp->myWeight = distance;
				dtp->myCount = label;
			}
		}
		label++;
	}
	return distanceMap;
}

Z3i::DigitalSet pointsToLinkWithJunction(const vector<WeightedPointCount<Z3i::Point>* >& distanceMapJunctions,
										 const Z3i::DigitalSet& endPoints,
										 int currentIndexJunctions) {
	typedef vector<WeightedPointCount<Z3i::Point>* >::value_type Scalar;

	Z3i::DigitalSet pointsToLink(endPoints.domain());
	for (Scalar wp : distanceMapJunctions) {
		int label = wp->myCount;
		if (label == currentIndexJunctions+1) {
			pointsToLink.insert(wp->myPoint);
		}
	}

	return pointsToLink;
}

vector<Z3i::DigitalSet> decompositionCCToSets(const Z3i::DigitalSet& aSet) {
	Z3i::Object26_6 obj(Z3i::dt26_6, aSet);
	vector<Z3i::Object26_6> cc;
	back_insert_iterator<vector<Z3i::Object26_6> > inserter(cc);
	obj.writeComponents(inserter);
	vector<Z3i::DigitalSet> ccSet;
	for (const auto& o : cc)
		ccSet.push_back(o.pointSet());
	return ccSet;
}


vector<Z3i::Point> thinCurves(const Z3i::DigitalSet& curve, const vector<Z3i::Point>& points) {
	vector<Z3i::Point> thinPoints;
	Z3i::Object26_6 obj(Z3i::dt26_6, curve);
	Z3i::Point previous;
	for (const Z3i::Point& p : points) {
		//If first point we add it (mixed connectivity in skeletal parts)
		if (p == *points.begin()) {
			thinPoints.push_back(p);
	    }
		else {
			vector<Z3i::Point> neighbors;
			back_insert_iterator<vector<Z3i::Point> > inserter(neighbors);
			obj.writeNeighbors(inserter, p);
			if (neighbors.size() < 2) {
				thinPoints.push_back(p);
			}
			else {
				//look up the half space for branching points
				Z3i::RealVector dirVector = (p - previous).getNormalized();
				trace.info() << dirVector << endl;
				DigitalPlane<Z3i::Space> digPlane(p, dirVector);
				bool add = true;
				for (const Z3i::Point& n : neighbors) {
					if (!digPlane.isPointAbove(n))
						add = false;
				}
				if (add)
					thinPoints.push_back(p);
				break;
			}
		}
		previous = p;
	}
	return thinPoints;
}

Z3i::DigitalSet curveThinning(const Z3i::DigitalSet& skeleton, const vector<vector<Z3i::Point> >& curves) {
	Z3i::DigitalSet endPoints(skeleton.domain());
	Z3i::DigitalSet curveSet(skeleton.domain());
	for (const vector<Z3i::Point> v : curves) {
		if (v.size() > 0) {
			Z3i::Point e1 = *v.begin();
			Z3i::Point e2 = *v.rbegin();
			endPoints.insert(e1);
			endPoints.insert(e2);

			curveSet.insert(v.begin(), v.end());
		}
	}

	Z3i::Object26_6 obj(Z3i::dt26_6, curveSet);

	int nb_simple = 0;
	do
    {
		Z3i::DigitalSet & S = obj.pointSet();
		std::queue<Z3i::DigitalSet::Iterator> Q;
		for ( Z3i::DigitalSet::Iterator it = S.begin(); it != S.end(); ++it )
			if ( obj.isSimple( *it ) && endPoints.find(*it) == endPoints.end() )
				Q.push( it );
		nb_simple = 0;
		while ( ! Q.empty() )
        {
			Z3i::DigitalSet::Iterator it = Q.front();
			Q.pop();
			if ( obj.isSimple( *it ) && endPoints.find(*it) == endPoints.end())
            {
				S.erase( *it );
				++nb_simple;
            }
        }
    }
	while ( nb_simple != 0 );

    return obj.pointSet();

}


void ensureSkeletonConnexity(Z3i::DigitalSet& skeletonPoints, const Z3i::DigitalSet& setVolume) {
	typedef Z3i::Object26_6 ObjectType;

	ObjectType objectVolume(Z3i::dt26_6, setVolume);
	ObjectType objectImage(Z3i::dt26_6, skeletonPoints);

	vector<ObjectType> skeletonCC;
	back_insert_iterator< std::vector<ObjectType> > inserterCC( skeletonCC );
	objectImage.writeComponents(inserterCC);
	bool objectsToLink = true;
	while (objectsToLink) {
		objectsToLink = false;
		double distance = numeric_limits<double>::max();
		Z3i::Point currentLink, otherLink;
		int indexDelete1, indexDelete2;
		for (size_t i = 0; i < skeletonCC.size(); i++) {
			Z3i::DigitalSet currentCC = skeletonCC[i].pointSet();
			for (size_t j = 0; j < skeletonCC.size(); j++) {
				if (i == j) continue;
				Z3i::DigitalSet otherCC = skeletonCC[j].pointSet();
				for (auto itCurrent = currentCC.begin(), itCurrentE = currentCC.end();
					 itCurrent != itCurrentE; ++itCurrent) {
					Z3i::Point current = *itCurrent;
					for (auto itOther = otherCC.begin(), itOtherE = otherCC.end(); itOther != itOtherE;
						 ++itOther) {
						Z3i::Point other = *itOther;
						double currentDistance = Z3i::l2Metric(current, other);
						if (
							currentDistance < distance) {
							distance = currentDistance;
							currentLink = current;
							otherLink = other;
							indexDelete1 = i;
							indexDelete2 = j;
							objectsToLink = true;
						}
					}
				}
			}
		}
		if (objectsToLink) {
			vector<Z3i::Point> link = SurfaceTraversal::AStarAlgorithm(objectVolume, currentLink, otherLink);
			Z3i::DigitalSet toKeep = skeletonCC[indexDelete1].pointSet();
			Z3i::DigitalSet toDelete = skeletonCC[indexDelete2].pointSet();
			toKeep.insert(link.begin(), link.end());
			toKeep.insert(toDelete.begin(),
						  toDelete.end());
			ObjectType toAdd(Z3i::dt26_6, toKeep);
			if (indexDelete1 == skeletonCC.size() - 1) {
				skeletonCC.pop_back();
				skeletonCC[indexDelete2] = skeletonCC[indexDelete1-1];
				skeletonCC.pop_back();
			}
			else if (indexDelete1 == skeletonCC.size() - 2) {
				skeletonCC[indexDelete2] = skeletonCC[indexDelete1+1];
				skeletonCC.pop_back();
				skeletonCC.pop_back();
			}
			else if (indexDelete2 == skeletonCC.size() - 1) {
				skeletonCC.pop_back();
				skeletonCC[indexDelete1] = skeletonCC[indexDelete2-1];
				skeletonCC.pop_back();
			}
		    else if (indexDelete2 == skeletonCC.size() - 2) {
				skeletonCC[indexDelete1] = skeletonCC[indexDelete2+1];
				skeletonCC.pop_back();
				skeletonCC.pop_back();
			}
			else {
				skeletonCC[indexDelete1] = skeletonCC[skeletonCC.size() - 1];
				skeletonCC[indexDelete2] = skeletonCC[skeletonCC.size() - 2];
				skeletonCC.pop_back();
				skeletonCC.pop_back();
			}
			skeletonCC.push_back(toAdd);
			skeletonPoints.insert(link.begin(), link.end());
		}
	}
	trace.info() << skeletonCC.size() << endl;
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
	typedef ImageContainerBySTLVector<Domain, double> ImageDouble;
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

	typedef DigitalPlane<Z3i::Space> DigitalPlane;
	typedef DigitalPlaneSet<Z3i::Space> DigitalPlaneSet;

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
		("skeleton,s", po::bool_switch()->default_value(false), "vol file (medial axis)")
		("delta,d", po::value<double>()->default_value(2), "delta for ball radius")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		("radiusInside,R", po::value<double>()->default_value(10), "radius of the ball inside voronoi cell")
		("radiusNeighbour,r", po::value<double>()->default_value(10), "radius of the ball for the neighbourhood")
		("thresholdFeature,T", po::value<double>()->default_value(0.1), "feature threshold")
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
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	double R = vm["radiusInside"].as<double>();
	double r = vm["radiusNeighbour"].as<double>();
	double delta = vm["delta"].as<double>();
	double thresholdFeature = vm["thresholdFeature"].as<double>();
	bool isDT = vm["skeleton"].as<bool>();

	QApplication application(argc,argv);
	Viewer3D<> viewer;
	viewer.show();

	Image volume = VolReader<Image>::importVol(inputFilename);
	Z3i::Domain domainVolume = volume.domain();
	Z3i::DigitalSet setVolume(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
												  thresholdMin-1, thresholdMax);

	Z3i::Object26_6 graph(Z3i::dt26_6, setVolume);

	trace.info() << setVolume.size() << endl;
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


	Metric l2;
	//SaddleComputer<DTL2, BackgroundPredicate> saddleComputer(setVolume, dt, backgroundPredicate, R, r, delta);
	//vector<pair<Z3i::Point, double> > pointToEigenValue = saddleComputer.mapPointToEigenvalue();
	// Z3i::DigitalSet branchingPoints = saddleComputer.extractSaddlePoints(setVolume);
	// vector<Z3i::Object26_6> objSaddle = saddleComputer.saddleConnectedComponents(branchingPoints);
 	// Z3i::DigitalSet maxCurvaturePoints = saddleComputer.saddlePointsToOnePoint<Matrix>(objSaddle);
	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin-1, thresholdMax);
	// Z3i::DigitalSet saddlePoints(domainVolume);
	// for (const pair<Z3i::Point, double>& pair : pointToEigenValue) {
	// 	if (pair.second > thresholdFeature)
	// 		saddlePoints.insert(pair.first);
	// }


	WeightedPointCount* currentPoint = *setVolumeWeighted.begin();
	double distanceMax = currentPoint->myWeight+delta;

	VCM vcm( R, ceil( r ), l2, true );
	vcm.init( setVolume.begin(), setVolume.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, r );


	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );

	int i = 0;
	int numberLeft = setVolumeWeighted.size();

	Z3i::RealPoint normal(0,0,1);
	Z3i::RealPoint previousNormal=normal, seedNormal = normal;

	Z3i::DigitalSet connectedComponent3D(domainVolume), seedConnectedComponent3D(domainVolume);
	Z3i::DigitalSet endPoints(domainVolume);
	Z3i::DigitalSet branchingParts(domainVolume);
	Z3i::RealPoint realCenter;
	Z3i::Point centerOfMass;
	Z3i::Point previousCenter, seedCenter;
	Z3i::DigitalSet branches(setVolume.domain());
	bool isNewSeed = true;
	vector<Z3i::DigitalSet> planes;
	vector<DigitalPlaneSet> shellPoints;

	map<Z3i::Point, int> junctionPoints;

	map<Z3i::Point, Z3i::RealPoint> pointToNormal;

	trace.beginBlock("Computing skeleton");
	//Main loop to compute skeleton (stop when no vol points left to process)
	while (numberLeft > 0)
	{
		if (numberLeft < 0 * setVolumeWeighted.size()) break;
		i++;
		trace.progressBar((setVolumeWeighted.size() - numberLeft), setVolumeWeighted.size());
		currentPoint->myProcessed = true;
		double radius = r;
		Point closestPointToCurrent = *min_element(vPoints.begin(), vPoints.end(), [&](const Point& one, const Point& two) {
					return Z3i::l2Metric(one, currentPoint->myPoint) < Z3i::l2Metric(two, currentPoint->myPoint);
			});
		//Distance transform value for VCM radius
		if (isDT) {
			if (dt(closestPointToCurrent) > dt(currentPoint->myPoint))
				radius = dt(closestPointToCurrent);
			else
				radius = dt(currentPoint->myPoint);
			radius += delta;
			if (radius > 0) {
				chi = KernelFunction(1.0, radius);
			}
 		}

		// Compute discrete plane
		connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
		 													 currentPoint->myPoint, normal,
		 													 0, radius, distanceMax, 26, true);


	    realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);

		Z3i::DigitalSet toMark(connectedComponent3D.domain());
		for (const Z3i::Point& p : connectedComponent3D) {
			if (Z3i::l2Metric(p, currentPoint->myPoint) <= radius)
				toMark.insert(p);
		}
		VCMUtil::markConnectedComponent3D(setVolumeWeighted, toMark, 0);
		for (const Z3i::Point& p : toMark) {
			pointToNormal[p] = normal;
		}
		//Center of mass computation
		if (realCenter != Z3i::RealPoint()) {
			centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
			pointToNormal[centerOfMass] = normal;
			double radiusCurrent = computeRadiusFromIntersection(volume, centerOfMass, normal, radius*6)+sqrt(3);
			double radiusCurrentMinus = computeRadiusFromIntersection(volume, centerOfMass, -normal, radius*6)+sqrt(3);

			double radiusShell = std::max(radiusCurrent, radiusCurrentMinus);

			if (Z3i::l2Metric(currentPoint->myPoint, realCenter) <= sqrt(3)
				) {
				Z3i::DigitalSet startingPoint(domainVolume);
				startingPoint.insert(centerOfMass);
				Z3i::DigitalSet shell = computeHalfShell(centerOfMass, Z3i::RealVector(0,0,0), startingPoint, setVolume, radiusShell);
				unsigned int nbCC = computeDegree(shell);
				int r = rand() % 256, b = rand() % 256, g = rand() % 256;
				if (nbCC >= 3) {
					Z3i::DigitalSet halfShellNormal1 = computeHalfShell(centerOfMass, normal, startingPoint, setVolume, radiusShell);
					Z3i::DigitalSet halfShellNormal2 = computeHalfShell(centerOfMass, -normal, startingPoint, setVolume, radiusShell);
					unsigned int nbCC1 = computeDegree(halfShellNormal1);
					unsigned int nbCC2 = computeDegree(halfShellNormal2);
					Ball<Z3i::Point> ball(centerOfMass, dt(centerOfMass));
					std::vector<Z3i::Point> points;
					DigitalPlane emptyPlane;
					DigitalPlaneSet digPlaneSet(emptyPlane, Z3i::DigitalSet(domainVolume));
					if (nbCC1 > nbCC2 && nbCC1 > 1 && nbCC2 > 0) {
						points = ball.pointsInHalfBall(normal);
						DigitalPlane digPlane(centerOfMass, -normal);
						digPlaneSet = DigitalPlaneSet(digPlane, halfShellNormal1);
						shellPoints.push_back(digPlaneSet);
					} else if (nbCC2 > nbCC1 && nbCC2 > 1 && nbCC2 > 0) {
						points = ball.pointsInHalfBall(-normal);
						DigitalPlane digPlane(centerOfMass, normal);
						digPlaneSet = DigitalPlaneSet(digPlane, halfShellNormal2);
						shellPoints.push_back(digPlaneSet);
					}
					junctionPoints[centerOfMass] = nbCC;

				}
				if (nbCC == 1) {
					// if ((nbCC1 == 0 && nbCC2 == 1) ||
					// 	(nbCC1 == 1 && nbCC2 == 0)) {
						endPoints.insert(centerOfMass);
					// }
				}
				if (Z3i::l2Metric(currentPoint->myPoint, realCenter) <= sqrt(3)) {
					bool processed = false;
					for (auto it = connectedComponent3D.begin(), ite = connectedComponent3D.end(); it != ite; ++it) {
						for (auto itS = skeletonPoints.begin(), itSe = skeletonPoints.end(); itS != itSe; ++itS)  {
							if (*itS == *it)
								processed = true;
						}
					}

					if (!processed && Z3i::l2Metric(currentPoint->myPoint, realCenter) <= sqrt(3)
						){
						if (!isNewSeed && Z3i::l2Metric(previousCenter, centerOfMass) <= 2 * sqrt(3)
							) {
							Z3i::DigitalSet diffPlanes = VCMUtil::markDifferenceBetweenPlanes(setVolumeWeighted,
																							  previousNormal, previousCenter,
																							  normal, centerOfMass,
																							  domainVolume, radius);
							planes.push_back(diffPlanes);
							VCMUtil::markConnectedComponent3D(setVolumeWeighted, diffPlanes, 0);

						}

						// Branching detection

						Z3i::Object26_6 objSkel(Z3i::dt26_6, skeletonPoints);
						vector<Z3i::Point> neighbors;
						back_insert_iterator<vector<Z3i::Point> > inserter(neighbors);
						objSkel.writeNeighbors(inserter, centerOfMass);
						if (neighbors.size() < 2) {
							skeletonPoints.insert(centerOfMass);
						}
						previousNormal = normal;
						previousCenter = centerOfMass;

						// viewer << CustomColors3D(Color::Red, Color::Red) << centerOfMass;
						// viewer << Viewer3D<>::updateDisplay;
						// qApp->processEvents();

					}
				}
			}
		}

			//Go to next point according to normal OR to max value in DT
			// if (isNewSeed) {
			// 	seedNormal = normal;
			// 	seedCenter = centerOfMass;
			// 	seedConnectedComponent3D = connectedComponent3D;
			// }
		bool previousSeed = isNewSeed;
		isNewSeed = VCMUtil::trackNextPoint(currentPoint, setVolumeWeighted, connectedComponent3D, centerOfMass, normal);
		if (!isNewSeed && !previousSeed)
			planes.push_back(connectedComponent3D);

		// if (isNewSeed) {
		// 	trace.info() << currentPoint->myPoint << endl;
		// 	WeightedPointCount* otherPoint = new WeightedPointCount(*currentPoint);
		// 	isNewSeed = VCMUtil::trackNextPoint(otherPoint, setVolumeWeighted, seedConnectedComponent3D, seedCenter, seedNormal);
		// 	if (!isNewSeed) {
		// 		currentPoint = otherPoint;
		// 	}
		// }
		numberLeft = count_if(setVolumeWeighted.begin(), setVolumeWeighted.end(),
							  [&](WeightedPointCount* wpc) {
								  return (!wpc->myProcessed);
							  });
	}
	trace.endBlock();

	fillHoles(skeletonPoints, dt);

	branches = shellPointsToShellAreas (setVolume, junctionPoints, skeletonPoints);
	Z3i::DigitalSet newBranches = dilateJunctions(setVolume, skeletonPoints, branches);
	for (const Z3i::Point& p : newBranches)
		junctionPoints[p] = 3;
	branches = shellPointsToShellAreas (setVolume, junctionPoints, skeletonPoints);
	//Discarding points being in branching parts
	for (auto it = branches.begin(), ite = branches.end(); it != ite; ++it) {
	  	auto itToErase = skeletonPoints.find(*it);
	  	if (itToErase != skeletonPoints.end())
	 		skeletonPoints.erase(itToErase);
	}

	Z3i::DigitalSet endPointsParts = endPointCurves(skeletonPoints, endPoints);
	vector<Z3i::DigitalSet> junctionCCSet = decompositionCCToSets(branches);
	vector<Z3i::DigitalSet> parts = decompositionCCToSets(skeletonPoints);
	orientNormalsEndPoints(pointToNormal, endPointsParts, parts);
	vector<WeightedPointCount*> distanceMapJunctions = distanceMapFromBranches (skeletonPoints, junctionCCSet, endPointsParts,
																				parts, pointToNormal, viewer);
	int maxLabel = (*max_element(distanceMapJunctions.begin(), distanceMapJunctions.end(),
								[&](WeightedPointCount* one, WeightedPointCount* two) {
									return one->myCount < two->myCount;
								 }))->myCount;

	trace.info() << junctionCCSet.size() << endl;

	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << skeletonPoints;
	for (int i = 0; i < junctionCCSet.size(); i++) {
		Z3i::DigitalSet pointsToLink = pointsToLinkWithJunction(distanceMapJunctions, endPointsParts, i);
		if (pointsToLink.size() == 0) continue;
		Z3i::Point ref = orientationNormal( pointsToLink, pointToNormal );
		vector<vector<Z3i::Point> > bezierCurves;
		for (const Z3i::Point& p : pointsToLink) {
			if (p == ref) continue;
			vector<Z3i::Point> curves = linkWithBezierCurves (pointToNormal, skeletonPoints, setVolumeWeighted, setVolume, ref, p, dt);
			bool add = true;
			for (const Z3i::Point& c : curves) {
				if (dt(c) == 0)
					add = false;
			}
// 			vector<Z3i::Point> thinPoints = thinCurves(skeletonPoints, curves);
// 			vector<Z3i::Point> thinPoints = curves;
			if (add) {
				bezierCurves.push_back(curves);
			}
		}
		Z3i::DigitalSet thinBezierCurves = curveThinning( skeletonPoints, bezierCurves);
		for (const Z3i::Point& c : thinBezierCurves) {
			skeletonPoints.insert(c);
		}
		int label = i;
		unsigned char r = (label * 80) % 256;
		unsigned char g = (label * 100) % 256;
		unsigned char b = (label * 50) % 256;
		viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << pointsToLink;
	}


	viewer << CustomColors3D(Color::Blue, Color::Blue) << skeletonPoints;
	ensureSkeletonConnexity(skeletonPoints, setVolume);
	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << skeletonPoints;
	for (auto it = branches.begin(), ite = branches.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color(0,50,0,50), Color(0,50,0,50)) <<*it;
	}

	//second pass
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
		(*it)->myProcessed = false;
	}

	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
	  	viewer << CustomColors3D(Color(0,0,120,10), Color(0,0,120,10)) << (*it)->myPoint;
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
