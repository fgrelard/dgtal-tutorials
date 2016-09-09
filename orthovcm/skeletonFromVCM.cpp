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
#include "DGtal/io/writers/ITKWriter.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/math/Histogram.h"
#include "DGtal/math/Statistic.h"
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
#include "geometry/SaddleComputer.h"
#include "clustering/diana.hpp"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "geometry/WeightedPointCount.h"
#include "surface/SurfaceTraversal.h"
#include "Statistics.h"
#include "shapes/Ball.h"
///////////////////////////////////////////////////////////////////////////////
#include "ReverseHatPointFunction.h"
#include "clustering/Watershed.h"

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
	Ball<Z3i::Point> ballOuterBall(center, radiusOuterBall);

	Z3i::DigitalSet shell(setVolume.domain());
	for (auto it = setVolume.begin(), ite = setVolume.end();
		 it != ite; ++it) {
		Z3i::Point current = *it;
		if (!(ballInnerBall.contains(current)) && ballOuterBall.contains(current)) {
			shell.insert(current);
		}
	}

	return shell;
}


unsigned int computeDegree(const Z3i::DigitalSet& shell) {
	typedef Z3i::Object26_6 ObjectType;

	ObjectType objectImage(Z3i::dt26_6, shell);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);

	return nbConnectedComponents;
}


template <typename VCM, typename Domain>
Z3i::DigitalSet computeBranchingPartsWithVCMFeature(const VCM& vcm,
													const Domain& domain, double threshold) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	Z3i::DigitalSet aSet(domain);
	double R = vcm.R();
	for (P2EConstIterator  it = vcm.mapPoint2ChiVCM().begin(),
			  itE = vcm.mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
		if (ratio > threshold &&  sqrt(lambda[2]) < R*R)
			aSet.insert(it->first);
	}
	return aSet;
}

template <typename VCM>
vector<double> extractCurvatureOnPoints(const VCM& vcm) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

    vector<double> curvatureValues;
	for (P2EConstIterator  it = vcm.mapPoint2ChiVCM().begin(),
			  itE = vcm.mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
	    curvatureValues.push_back(ratio);
	}
	return curvatureValues;
}

template <typename DTL2>
Z3i::Point projectionOnSurface(const DTL2& dt, const Z3i::DigitalSet& aSet, const Z3i::Point& origin, const Z3i::RealPoint& normal) {
	int scalar = 1;
	Z3i::Point current = origin + normal;
	while (aSet.find(current) != aSet.end() && (current != origin || dt(current) > 1)) {
		scalar++;
		current = origin + normal * scalar;
	}
	return current;
}

template <typename VCM>
Z3i::RealPoint extractVectorVCMAtPoint(const VCM& vcm, const Z3i::Point& point, int coordinate) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	auto mapPointEigen = vcm.mapPoint2ChiVCM();
	Z3i::RealPoint vector = mapPointEigen[point].vectors.column(coordinate);
	double scalar = mapPointEigen[point].values[coordinate];
	return vector;
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
	Z2i::DigitalSet aSet = VCMUtil::extractConnectedComponent(processImage, center, 1, 255);
	Eigen::MatrixXd covmatrix = Statistics::computeCovarianceMatrix<Eigen::MatrixXd>(aSet);
	if (covmatrix.size() == 0) return 0;
	Z2i::RealVector projection = Statistics::extractEigenVector<Z2i::RealVector>(covmatrix, 1);
	Z2i::Point trackedPoint = PointUtil::trackPoint(center, aSet, projection);
	double distance = Z2i::l2Metric(center, trackedPoint);
	return distance;
	// if (covmatrix.size() == 0) return 0;
	// double eigenvalue = Statistics::extractEigenValue<Z2i::RealPoint>(covmatrix, 0)[1];
	// if (eigenvalue >= 0)
	// 	return sqrt(eigenvalue);
	// else
	// 	return 0;
}

template <typename Matrix>
double computeRadiusFromIntersection3D(const Z3i::DigitalSet& aSet, const Z3i::Point& centerOfMass) {
	Matrix covmatrix = Statistics::computeCovarianceMatrix<Matrix>(aSet);
	if (covmatrix.size() == 0) return 0;
	double eigenvalue = Statistics::extractEigenValue<Z3i::RealPoint>(covmatrix, 0)[2];
	if (eigenvalue >= 0)
		return sqrt(5.991*eigenvalue);
	else
		return 0;
}




Z3i::DigitalSet detectBranchingPointsInNeighborhood(Z3i::Point& branchingPoint,
													const Z3i::DigitalSet& branchingPoints,
													const Z3i::Point& current, double radius) {

	Z3i::DigitalSet aSet(branchingPoints.domain());
	Ball<Z3i::Point> ball(current, radius);

//	Z3i::Point branchingPoint;
	double distanceMin = numeric_limits<double>::max();
	bool toMark = false;
	for (auto it = branchingPoints.begin(), ite = branchingPoints.end(); it != ite; ++it) {
		if (ball.contains(*it) &&
			Z3i::l2Metric(*it, current) < distanceMin) {
			toMark = true;
			branchingPoint = *it;
		    distanceMin = Z3i::l2Metric(*it, current);
		}
	}
	if (toMark) {
		double ballRadius = Z3i::l2Metric(branchingPoint, current) + 1;
//		ballRadius = radius;
		Ball<Z3i::Point> ballBranching(current, ballRadius);
		Z3i::RealPoint dirVector = ((Z3i::RealPoint) branchingPoint - (Z3i::RealPoint)current).getNormalized();
		std::vector<Z3i::Point> pointsInBall;
		if (branchingPoint != current)
			 pointsInBall = ballBranching.pointsInHalfBall(dirVector);
		else
			pointsInBall = ballBranching.pointsInBall();
		aSet.insert(pointsInBall.begin(), pointsInBall.end());
	}
	return aSet;
}



template <typename VCM>
double radiusForJunction(const Z3i::DigitalSet& branchingPoints, const VCM& vcm,
													const Z3i::RealPoint& current, double radius) {
	Z3i::Point nearestPointBranch;
	double minDistance = numeric_limits<double>::max();
	Z3i::DigitalSet aSet(branchingPoints.domain());
	for (auto it = branchingPoints.begin(), ite = branchingPoints.end(); it != ite; ++it) {
		if (Z3i::l2Metric(*it, current) <= minDistance) {
			nearestPointBranch = *it;
			minDistance = Z3i::l2Metric(*it, current);
		}
	}
	if (nearestPointBranch != Z3i::Point()) {
		return VCMUtil::radiusAtJunction(vcm, nearestPointBranch, radius);
	}
	return 0;

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


Z3i::DigitalSet neighborsOfSet(const Z3i::DigitalSet& aSet, const Z3i::DigitalSet& constraintSet) {

	Z3i::DigitalSet neighbors(aSet.domain());
	DigitalSetInserter<Z3i::DigitalSet> inserter(neighbors);
	for (auto it = aSet.begin(), ite = aSet.end(); it != ite; ++it) {
		Z3i::Adj26::writeNeighbors(inserter, *it);
	}
	Z3i::DigitalSet neighborsInConstraintSet(aSet.domain());
	for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
		if (constraintSet.find(*it) != constraintSet.end())
			neighborsInConstraintSet.insert(*it);
	}
	return neighborsInConstraintSet;
}





template <typename VCM, typename KernelFunction>
double computeCurvature(const Z3i::DigitalSet& junctions, const Z3i::Point& branchPoint, double radius=1.0) {
	VCM vcm( radius, ceil( radius ), Z3i::l2Metric, false );
	vcm.init( junctions.begin(), junctions.end() );
	KernelFunction chi( 1.0, radius );
	Z3i::RealPoint lambda = VCMUtil::computeEigenValuesFromVCM(branchPoint, vcm, chi);
	return VCMUtil::computeCurvatureJunction(lambda);
}



Z3i::Point nearestBranchPoint(const Z3i::DigitalSet& aSet, const Z3i::DigitalSet& branchingPoints) {
	Z3i::Point point = *(aSet.begin());
	Z3i::Point closestPoint = *min_element(branchingPoints.begin(), branchingPoints.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
			return Z3i::l2Metric(one, point) < Z3i::l2Metric(two, point);
		});
	return closestPoint;
}




template <typename Viewer>
void fillHoles(Z3i::DigitalSet& skeletonPoints, Viewer& viewer) {
	Z3i::Object26_6 objectImage(Z3i::dt26_6, skeletonPoints);
	vector<Z3i::Object26_6> skeletonCC;
	back_insert_iterator< std::vector<Z3i::Object26_6> > inserterCC( skeletonCC );
	objectImage.writeComponents(inserterCC);
	for (size_t i = 0; i < skeletonCC.size(); i++) {
		Z3i::DigitalSet currentCC = skeletonCC[i].pointSet();
		for (size_t j = i+1; j < skeletonCC.size(); j++) {
			Z3i::DigitalSet otherCC = skeletonCC[j].pointSet();
			double distance = numeric_limits<double>::max();
			Z3i::Point currentLink, otherLink;
			for (auto itCurrent = currentCC.begin(), itCurrentE = currentCC.end();
				 itCurrent != itCurrentE; ++itCurrent) {
				Z3i::Point current = *itCurrent;
				for (auto itOther = otherCC.begin(), itOtherE = otherCC.end(); itOther != itOtherE;
					 ++itOther) {
					Z3i::Point other = *itOther;
					double currentDistance = Z3i::l2Metric(current, other);
					if (currentDistance > sqrt(3) &&
						currentDistance <= 2* sqrt(3) &&
						currentDistance < distance) {
						distance = currentDistance;
						currentLink = current;
						otherLink = other;
					}
				}
			}
			vector<Z3i::Point> link = PointUtil::linkTwoPoints(currentLink, otherLink);
			for (auto itL = link.begin(), itLe = link.end(); itL != itLe; ++itL) {
				viewer << CustomColors3D(Color::Yellow, Color::Yellow) << *itL;
				skeletonPoints.insert(*itL);
			}
		}
	}
}



template <typename DT>
Z3i::DigitalSet erodeSaddleAreas(const Z3i::DigitalSet& saddles,
					  const Z3i::DigitalSet& surface,
					  const vector<Z3i::Point>& medialAxis,
					  const DT& dt) {
	Z3i::DigitalSet saddleEroded(saddles.domain());
	for (const Z3i::Point& saddlePoint : saddles) {
		Z3i::Point closestPointToCurrent = *min_element(medialAxis.begin(), medialAxis.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
				return Z3i::l2Metric(one, saddlePoint) < Z3i::l2Metric(two, saddlePoint);
			});
		double radius = dt(closestPointToCurrent);
		Ball<Z3i::Point> ball(saddlePoint, 1);
		vector<Z3i::Point> surfaceIntersection = ball.surfaceIntersection(surface);
		bool add = true;
		for (const Z3i::Point& pointSI : surfaceIntersection) {
			if (saddles.find(pointSI) == saddles.end()) {
				add=false;
			}
		}
		if (add)
			saddleEroded.insert(saddlePoint);
	}
	return saddleEroded;
}


map<Z3i::Point, set<Z3i::Point> > shortestPathsJunctions(const Z3i::DigitalSet& setTraversed,
									   const Z3i::DigitalSet& setOfPoints,
									   const Z3i::Point& referencePoint) {
	Z3i::Point refPrime = extractNearestNeighborInSetFromPoint(setTraversed, referencePoint);
	Z3i::Object26_6 obj(Z3i::dt26_6, setTraversed);
	map<Z3i::Point, set<Z3i::Point>> setPaths;
	for (const Z3i::Point& p : setOfPoints) {
		vector<Z3i::Point> path = SurfaceTraversal::AStarAlgorithm(obj, p, refPrime);
		setPaths[p].insert(path.begin(), path.end());
	}
	return setPaths;
}

Z3i::DigitalSet ensureConnexity(const Z3i::DigitalSet& set) {
	Z3i::DigitalSet cleanSet(set.domain());
	Z3i::Object26_6 obj(Z3i::dt26_6, set);
	Z3i::DigitalSet & S = obj.pointSet();
	cleanSet = S;
	for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
		Z3i::Object26_6 objClean(Z3i::dt26_6, cleanSet);
		if (objClean.isSimple(*it)) {
		    cleanSet.erase(*it);
		}
	}

	return cleanSet;
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



template <typename VCM, typename KernelFunction>
Z3i::DigitalSet jordanCurve(const VCM& vcm, const KernelFunction& chi,
							const map<Z3i::Point, set<Z3i::Point> >& paths,
							const Z3i::DigitalSet& setTraversed,
							const Z3i::DigitalSet& setOfPoints) {

	typedef Eigen::MatrixXd Matrix;
	Z3i::Object26_6 obj(Z3i::dt26_6, setOfPoints);
	Z3i::Object26_6 objSurface(Z3i::dt26_6, setTraversed);
	Z3i::DigitalSet jordan(setTraversed.domain());
	set<Z3i::Point> setSurface(setTraversed.begin(), setTraversed.end());
	int nb = 0;
	set<Z3i::Point> processedSet;
	for (auto it = paths.begin(), ite = paths.end(); it != ite; ++it) {
		if (processedSet.find(it->first) != processedSet.end()) continue;
		Z3i::Point p = it->first;
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
		obj.writeNeighbors(inserter, p);
		set<Z3i::Point> pathCurrent = paths.at(p);
		set<Z3i::Point> unionPath;

		for (const Z3i::Point& n : neighbors) {
			if (processedSet.find(n) != processedSet.end()) {
				unionPath = pathCurrent;
				break;
			}
			set<Z3i::Point> pathNeighbor = paths.at(n);
			set_union(pathCurrent.begin(), pathCurrent.end(),
					  pathNeighbor.begin(), pathNeighbor.end(),
					  std::inserter(unionPath, unionPath.end()));

			Z3i::DigitalSet unionDigitalSet(setTraversed.domain());
			unionDigitalSet.insert(unionPath.begin(), unionPath.end());
			Z3i::Domain domain = PointUtil::computeBoundingBox<Z3i::Domain>(unionPath);
			Z3i::Point center = Statistics::extractCenterOfMass3D(unionDigitalSet);
			Matrix covMatrix = Statistics::computeCovarianceMatrix<Matrix>(unionDigitalSet);
			if (covMatrix.size() == 0) continue;
			Z3i::RealVector normalPlane = Statistics::extractEigenVector<Z3i::RealVector>(covMatrix, 0);
			set<Z3i::Point> intersectionPlaneHull;
			double omega = std::abs(normalPlane[0]) + std::abs(normalPlane[1]) + std::abs(normalPlane[2]);
			double d = -(-normalPlane[0] * center[0] - normalPlane[1] * center[1] - normalPlane[2] * center[2]);

			for (const Z3i::Point& p : domain) {
				double valueToCheckForPlane = p[0] * normalPlane[0] + p[1] * normalPlane[1] + p[2] * normalPlane[2];
				if (valueToCheckForPlane >= d && valueToCheckForPlane < d + omega)
					intersectionPlaneHull.insert(p);
			}

			set<Z3i::Point> cut;
			set_intersection(intersectionPlaneHull.begin(), intersectionPlaneHull.end(),
							 setSurface.begin(), setSurface.end(),
							 std::inserter(cut, cut.end()));
//		trace.info() << cut.size() << " " << normalPlane << intersectionPlaneHull.size() << endl;

			set<Z3i::Point> domainSet;
			domainSet.insert(domain.begin(), domain.end());
			set<Z3i::Point> difference;
			set_difference(setSurface.begin(), setSurface.end(),
						   domainSet.begin(), domainSet.end(),
						   std::inserter(difference, difference.end()));
			Z3i::DigitalSet diffDigitalSet(setTraversed.domain());
			diffDigitalSet.insert(difference.begin(), difference.end());
			Z3i::Object26_6 objDiff(Z3i::dt26_6, diffDigitalSet);
			vector<Z3i::Object26_6> cc;
			back_insert_iterator< vector<Z3i::Object26_6> > iterator( cc );
			unsigned nbcc = objDiff.writeComponents(iterator);
			if (nbcc >= 2
				) {
				processedSet.insert(neighbors.begin(), neighbors.end());
				processedSet.insert(p);
				jordan.insert(unionPath.begin(), unionPath.end());
				nb++;
			}
			unionPath.clear();
		}



		// set<Z3i::Point> upath = unionPath;
		// for (const Z3i::Point& up : unionPath) {
		// 	Ball<Z3i::Point> ball(up, 2);
		// 	vector<Z3i::Point> surfaceIntersection = ball.surfaceIntersection(setTraversed);
		// 	upath.insert(surfaceIntersection.begin(), surfaceIntersection.end());
		// }
		// Z3i::DigitalSet unionPathSet(setTraversed.domain());
		// unionPathSet.insert(unionPath.begin(), unionPath.end());


		// set<Z3i::Point> difference;
		// set_difference(setSurface.begin(), setSurface.end(),
		// 			   upath.begin(), upath.end(),
		// 			   std::inserter(difference, difference.end()));

		// Z3i::DigitalSet diffDigitalSet(setTraversed.domain());
		// diffDigitalSet.insert(difference.begin(), difference.end());
		// Z3i::Object26_6 objDiff(Z3i::dt26_6, diffDigitalSet);
		// vector<Z3i::Object26_6> cc;
		// back_insert_iterator< vector<Z3i::Object26_6> > iterator( cc );
		// unsigned nbcc = objDiff.writeComponents(iterator);
		// if (nbcc >= 2
		// 	) {
		// 	processedSet.insert(neighbors.begin(), neighbors.end());
		// 	processedSet.insert(p);
		// 	jordan.insert(upath.begin(), upath.end());
		// 	nb++;
		// }
	}
	trace.info() << nb << endl;
	if (nb >= 2)
		return jordan;
	else
		return Z3i::DigitalSet(setTraversed.domain());
}

bool planesIntersect(const Z3i::DigitalSet& plane,
					 const Z3i::DigitalSet& plane2) {
	for (const Z3i::Point& p : plane) {
		if (plane2.find(p) != plane2.end())
			return true;
	}
	return false;
}

vector<pair<Z3i::DigitalSet, vector<Z3i::DigitalSet> > > intersectingPlanes(const vector<Z3i::DigitalSet>& planes) {
	vector<pair<Z3i::DigitalSet, vector<Z3i::DigitalSet> > > inter;
	for (int i = 0, end = planes.size(); i < end; i++) {
		Z3i::DigitalSet current = planes[i];
		vector<Z3i::DigitalSet> planesIntersecting;
		for (int j = i+1; j < end; j++) {
			Z3i::DigitalSet other = planes[j];
			if (planesIntersect(current, other))
				planesIntersecting.push_back(other);
		}
		pair<Z3i::DigitalSet, vector<Z3i::DigitalSet>> pair = make_pair(current, planesIntersecting);
		inter.push_back(pair);
	}
	return inter;
}

double highestAreaVariation(const pair<Z3i::DigitalSet, vector<Z3i::DigitalSet>>& pairPlanes) {
	double currentArea = pairPlanes.first.size();
	double factorMax = 0.0;

	for (const Z3i::DigitalSet& otherPlane : pairPlanes.second) {
		double otherArea = otherPlane.size();
		double factor = otherArea / currentArea;
		if (factor > factorMax) {
			factorMax = factor;
		}
	}
	return factorMax;
}

template <typename VCM, typename KernelFunction, typename Image, typename Container>
pair<Z3i::DigitalSet, double> eigenValuesWithVCM(VCM& vcm, KernelFunction& chi, const Image& volume,
												 const Z3i::Point& p, double radius,
												 const Container& setVolumeWeighted, const Z3i::DigitalSet& setVolume,
												 Viewer3D<>& viewer) {
	map<Z3i::DigitalSet, double> eigenValues;
	Z3i::RealPoint normal;
	Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, setVolume.domain(), setVolumeWeighted,
																		 p, normal,
																		 0, radius, radius*10, 26, true);

	//  Z3i::RealPoint realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);
    //  Z3i::Point	centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
	//  double radiusIntersection = radius*2;
	// connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, setVolume.domain(), setVolumeWeighted,
	// 													 centerOfMass, normal,
	// 													 0, radiusIntersection, radius*10, 26, false);

	Z3i::RealVector eigenValue = VCMUtil::computeEigenValuesFromVCM(p, vcm, chi);
	vector<Z3i::Point> v(connectedComponent3D.begin(), connectedComponent3D.end());
	double eigenVar = sqrt(eigenValue[0]);
	double angle = M_PI/2 - vcm.vectorVariability(v, normal, chi, p);
	// Z3i::DigitalSet shell = computeShell (p, setVolume, radius*2, radius*4);
	// int degree = computeDegree(shell);
	// if (degree == 1)
	// 	angle = 0;
	Z3i::DigitalSet discretePlane(connectedComponent3D.domain());
	for (const Z3i::Point& c : connectedComponent3D) {
		if (Z3i::l2Metric(p, c) <= radius)
			discretePlane.insert(c);
	}

	return make_pair(discretePlane, angle);
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
	// SaddleComputer<DTL2, BackgroundPredicate> saddleComputer(setVolume, dt, backgroundPredicate, R, r, delta);
	// Z3i::DigitalSet branchingPoints = saddleComputer.extractSaddlePoints(setVolume);
	// vector<Z3i::Object26_6> objSaddle = saddleComputer.saddleConnectedComponents(branchingPoints);
 	// Z3i::DigitalSet maxCurvaturePoints = saddleComputer.saddlePointsToOnePoint<Matrix>(objSaddle);
	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin-1, thresholdMax);



	WeightedPointCount* currentPoint = *setVolumeWeighted.begin();
	double distanceMax = currentPoint->myWeight+delta;

	VCM vcm( R, ceil( r ), l2, true );
	viewer << CustomColors3D(Color::Green, Color::Green);
	vcm.init( setVolume.begin(), setVolume.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, r );


	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );

	int i = 0;
	int numberLeft = setVolumeWeighted.size();

	Z3i::RealPoint normal(0,0,1);
	Z3i::RealPoint previousNormal=normal, seedNormal = normal;

	Z3i::DigitalSet connectedComponent3D(domainVolume), seedConnectedComponent3D(domainVolume);
	Z3i::DigitalSet branchingParts(domainVolume);
	Z3i::RealPoint realCenter;
	Z3i::Point centerOfMass;
	Z3i::Point previousCenter, seedCenter;
	Z3i::DigitalSet branches(setVolume.domain());
	bool isNewSeed = true;
	vector<Z3i::DigitalSet> planes;
    map<Z3i::Point, double> pointToEigenValue;

	for (const Z3i::Point& p : setVolume) {
		pointToEigenValue[p] = numeric_limits<double>::max();
	}
	trace.beginBlock("Computing skeleton");

	//Main loop to compute skeleton (stop when no vol points left to process)
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it)
	//while (numberLeft > 0)
	{
		i++;
		trace.progressBar(i, setVolumeWeighted.size());
		//trace.progressBar((setVolumeWeighted.size() - numberLeft), setVolumeWeighted.size());
		currentPoint = *it;
		currentPoint->myProcessed = true;
		double radius = r;

		//Distance transform value for VCM radius
		if (isDT) {
			Point closestPointToCurrent = *min_element(vPoints.begin(), vPoints.end(), [&](const Point& one, const Point& two) {
					return Z3i::l2Metric(one, currentPoint->myPoint) < Z3i::l2Metric(two, currentPoint->myPoint);
				});
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
		// connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
		// 													 currentPoint->myPoint, normal,
		// 													 0, radius, distanceMax, 26, true);

	    pair<Z3i::DigitalSet, double> value = eigenValuesWithVCM (vcm, chi, volume, currentPoint->myPoint, radius, setVolumeWeighted, setVolume, viewer);
		pointToEigenValue[currentPoint->myPoint] =  std::min(pointToEigenValue[currentPoint->myPoint], value.second);
		connectedComponent3D = value.first;
		// for (const Z3i::Point& p : value.first) {
		//     if (connectedComponent3D.find(p) != connectedComponent3D.end())
		// 		pointToEigenValue[p] = std::min(pointToEigenValue[p], value.second);
		// }

	    // realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);


		// //Center of mass computation
		// if (realCenter != Z3i::RealPoint()) {
		// 	centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);

		// 	double radiusIntersection = computeRadiusFromIntersection(volume, centerOfMass, normal, distanceMax*2);

		// 	bool processed = false;
		// 	for (auto it = connectedComponent3D.begin(), ite = connectedComponent3D.end(); it != ite; ++it) {
		// 		for (auto itS = skeletonPoints.begin(), itSe = skeletonPoints.end(); itS != itSe; ++itS)  {
		// 			if (*itS == *it)
		// 				processed = true;
		// 		}
		// 	}

		// 	//VCMUtil::markConnectedComponent3D(setVolumeWeighted, connectedComponent3D, 0);

		// 	if (!processed && Z3i::l2Metric(currentPoint->myPoint, centerOfMass) <= sqrt(3)
		// 		){
		// 		if (!isNewSeed && Z3i::l2Metric(previousCenter, centerOfMass) <= 2 * sqrt(3)
		// 			) {
		// 			// Z3i::DigitalSet diffPlanes = VCMUtil::markDifferenceBetweenPlanes(setVolumeWeighted,
		// 			// 																  previousNormal, previousCenter,
		// 			// 																   normal, centerOfMass,
		// 			// 																   domainVolume, radius);
		// 			// planes.push_back(diffPlanes);
		// 			// for (const Z3i::Point & p : diffPlanes) {
		// 			// 	if (DGtal::Z3i::l2Metric(p, currentPoint->myPoint) <= radius)
		// 			// 		pointToEigenValue[p] = std::min(pointToEigenValue[p], value.second);
		// 			// }
		// 			//	VCMUtil::markConnectedComponent3D(setVolumeWeighted, diffPlanes, 0);

		// 		}

		// 		// Branching detection
		// 		skeletonPoints.insert(centerOfMass);
		// 		previousNormal = normal;
		// 		previousCenter = centerOfMass;

		// 		viewer << CustomColors3D(Color::Red, Color::Red) << centerOfMass;
		// 		viewer << Viewer3D<>::updateDisplay;
		// 		qApp->processEvents();

		// 	}
		// }

		// //Go to next point according to normal OR to max value in DT
		// // if (isNewSeed) {
		// // 	seedNormal = normal;
		// // 	seedCenter = centerOfMass;
		// // 	seedConnectedComponent3D = connectedComponent3D;
		// // }
		// bool previousSeed = isNewSeed;
		// isNewSeed = VCMUtil::trackNextPoint(currentPoint, setVolumeWeighted, connectedComponent3D, centerOfMass, normal);
		// if (!isNewSeed && !previousSeed)
		// 	planes.push_back(connectedComponent3D);

		// // if (isNewSeed) {
		// // 	trace.info() << currentPoint->myPoint << endl;
		// // 	WeightedPointCount* otherPoint = new WeightedPointCount(*currentPoint);
		// // 	isNewSeed = VCMUtil::trackNextPoint(otherPoint, setVolumeWeighted, seedConnectedComponent3D, seedCenter, seedNormal);
		// // 	if (!isNewSeed) {
		// // 		currentPoint = otherPoint;
		// // 	}
		// // }
		// numberLeft = count_if(setVolumeWeighted.begin(), setVolumeWeighted.end(),
		// 					  [&](WeightedPointCount* wpc) {
		// 						  return (!wpc->myProcessed);
		// 					  });
  		// i++;
	}
	trace.endBlock();
	double minVal = min_element(pointToEigenValue.begin(), pointToEigenValue.end(), [&](const pair<Z3i::Point, double>& one,
																						const pair<Z3i::Point, double>& two) {
									return (one.second < two.second);
								})->second;
	double maxVal = 0;
	for (const auto & pair : pointToEigenValue) {
		if (pair.second == numeric_limits<double>::max()) continue;
		if (maxVal < pair.second)
			maxVal = pair.second;
	}
	trace.info() << minVal << " " << maxVal << endl;
	vector<map<Z3i::Point, double>::iterator> itToRemove;
	for (auto it = pointToEigenValue.begin(), ite = pointToEigenValue.end(); it != ite; ++it) {
		if (it->second == numeric_limits<double>::max())
			itToRemove.push_back(it);
	}
	for (const auto& iterator : itToRemove) {
	    pointToEigenValue.erase(iterator);
	}
	// trace.info() << pointToEigenValue.size() << " " << setVolume.size() << endl;
	// trace.info() << minVal << " " << maxVal << endl;
	// Watershed<Z3i::Point> watershed(pointToEigenValue, thresholdFeature);
	// watershed.compute();
	// auto resultWatershed = watershed.getWatershed();
    // int bins = watershed.getBins();

	// GradientColorMap<int, CMAP_JET > hueShade(0, bins);
	// trace.info() << "Bins= " << bins << endl;
	// for (const auto& complexPoint : resultWatershed) {
	// 	viewer << CustomColors3D(hueShade(complexPoint.second.getLabel()), hueShade(complexPoint.second.getLabel())) << complexPoint.first;
	// }

	GradientColorMap<double, CMAP_JET > hueShade(minVal, maxVal);
	for (const auto& pToL : pointToEigenValue) {
		if (pToL.second > maxVal) continue;
		viewer << CustomColors3D(hueShade(pToL.second), hueShade(pToL.second)) << pToL.first;
	}



	//Discarding points being in branching parts
	for (auto it = branches.begin(), ite = branches.end(); it != ite; ++it) {
	 	auto itToErase = skeletonPoints.find(*it);
	 	if (itToErase != skeletonPoints.end())
	 		skeletonPoints.erase(itToErase);
	}

//	fillHoles(skeletonPoints, viewer);

	// for (auto it=skeletonPoints.begin(), ite = skeletonPoints.end(); it != ite; ++it) {
	// 	viewer << CustomColors3D(Color::Red, Color::Red) << *it;
	// }

	// for (auto it = branches.begin(), ite = branches.end(); it != ite; ++it) {
	//  	viewer << CustomColors3D(Color(0,50,0,50), Color(0,50,0,50)) <<*it;
	// }

	//second pass
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
		(*it)->myProcessed = false;
	}
	// Z3i::Object26_6 objectifyVolume(Z3i::dt26_6, setVolume);
	// Z3i::DigitalSet junctions = connectDisconnectedComponents<VCM, KernelFunction>(objectifyVolume, skeletonPoints, branches, *vcm_surface, branchingPoints);
	// for (auto it = junctions.begin(), ite = junctions.end(); it != ite; ++it) {
	//  	viewer << CustomColors3D(Color::Green, Color::Green) << *it;
	// }



	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
	  	viewer << CustomColors3D(Color(0,0,120,10), Color(0,0,120,10)) << (*it)->myPoint;
	}

	Image outImage(volume.domain());
	ImageDouble outDoubleImage(volume.domain());
	for (const Z3i::Point& p : domainVolume) {
		if (setVolume.find(p) == setVolume.end())
			outDoubleImage.setValue(p, minVal-1);
	}
	for (const auto& pToL : pointToEigenValue) {
	    if (pToL.second > maxVal) continue;
		outDoubleImage.setValue(pToL.first, pToL.second);
	}
	DGtal::imageFromRangeAndValue(skeletonPoints.begin(), skeletonPoints.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);
	typedef ITKWriter<ImageDouble> Writer;
	Writer::exportITK(outFilename + ".mhd", outDoubleImage);
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
