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
//#include "surface/SurfaceUtils.h"
#include "surface/Morphomaths.h"
#include "clustering/diana.hpp"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "geometry/WeightedPointCount.h"
#include "surface/SurfaceTraversal.h"
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
	Ball<Z3i::Point> ballOuterBall(center, radiusOuterBall);

	Z3i::DigitalSet pointsInOuterBall = ballOuterBall.pointsInBallSet();

	Z3i::DigitalSet shell(pointsInOuterBall.domain());
	for (auto it = pointsInOuterBall.begin(), ite = pointsInOuterBall.end();
		 it != ite; ++it) {
		if (!(ballInnerBall.contains(*it)) && setVolume.find(*it) != setVolume.end()) {
			shell.insert(*it);
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
Z3i::DigitalSet computeBranchingPartsWithVCMFeature(const VCM& vcm, const Domain& domain, double threshold) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	Z3i::DigitalSet aSet(domain);
	for (P2EConstIterator  it = vcm.mapPoint2ChiVCM().begin(),
			  itE = vcm.mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
		if (ratio > threshold)
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

template <typename Container>
double otsuThreshold(const Container& container) {

	double  proba = 0;                // first order cumulative
    double  mu = 0;                // second order cumulative
    double  mean = 0;               // total mean level        
	double    threshold = 0;        // optimal threshold value
	double max = 0.0;
	
	Statistic<double> stats;
	stats.addValues( container.begin(), container.end() );
	stats.terminate(); // stats are computed.
	Histogram<double>* hist = new Histogram<double>();
	hist->init( Histogram<double>::SquareRoot, stats );
	hist->addValues( container.begin(), container.end() );
	hist->terminate();
	
	for (int i=0; i< hist->size(); i++) {
		mean+= ((double) i / hist->size()) * hist->pdf(i);
	}
	for (int i = 0; i < hist->size(); i++) {
		proba += hist->pdf(i);
		mu += ((double)i/hist->size()) * hist->pdf(i);
		double currentValue =  pow((mean * proba - mu), 2) * proba * (1 - proba);
		if (currentValue > max) {
			max = currentValue;
			threshold = ((double)i/hist->size());
		}
			
	}
		
	return threshold;
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

Z3i::DigitalSet markAroundBranchingPoints(const Z3i::DigitalSet& setVolume, const Z3i::Point& center, double radius) {
	Z3i::DigitalSet aSet(setVolume.domain());
	Ball<Z3i::Point> ball(center, radius);
	std::vector<Z3i::Point> pointsInBall = ball.pointsInBall();
	for (auto it = pointsInBall.begin(), ite = pointsInBall.end(); it != ite; ++it) {
		if (setVolume.find(*it) != setVolume.end()) {
			aSet.insert(*it);
		}
	}
	return aSet;
}

template <typename ImageAdapterExtractor, typename Matrix, typename Image>
double computeRadiusFromIntersection(const Image& volume, const Z3i::Point& point, const Z3i::RealPoint& normal,
									 double radius) {
	typedef ImageSelector<Z2i::Domain, bool>::Type Image2D;
	DGtal::functors::Identity idV;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-radius, -radius, -radius), volume.domain().upperBound() + Z3i::Point(radius, radius, radius));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
									  DGtal::Z2i::Point(radius, radius));
    
	
	DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, point, normal, radius, domain3Dyup.lowerBound());

	ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
	Image2D processImage = ImageUtil::convertImage<Image2D>(extractedImage);
	Z2i::DigitalSet aSet = VCMUtil::extractConnectedComponent(processImage, Z2i::Point(radius/2, radius/2), 1, 255);
    Matrix covmatrix = Statistics::computeCovarianceMatrix<Matrix>(aSet);
	if (covmatrix.size() == 0) return 0;
	double eigenvalue = Statistics::extractEigenValue<Z2i::RealPoint>(covmatrix, 0)[1];
	if (eigenvalue >= 0)
		return sqrt(eigenvalue);
	else
		return 0;
}

template <typename Matrix>
double computeRadiusFromIntersection3D(const Z3i::DigitalSet& aSet) {
	Matrix covmatrix = Statistics::computeCovarianceMatrix<Matrix>(aSet);
	if (covmatrix.size() == 0) return 0;
	double eigenvalue = Statistics::extractEigenValue<Z3i::RealPoint>(covmatrix, 0)[2];
	if (eigenvalue >= 0)
		return sqrt(eigenvalue);
	else
		return 0;
}

template <typename WeightedPoint>
set<WeightedPoint*> computePointsBelowPlane(set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume, const Z3i::RealPoint& normal, const Z3i::Point& center, double distanceMax = numeric_limits<double>::max()) {
	set<WeightedPoint*> inferiorPoints;
	double d = -(-normal[0] * center[0] - normal[1] * center[1] - normal[2] * center[2]);
	for (auto it = volume.begin(), ite = volume.end(); it!=ite; ++it) {
		double valueToCheckForPlane = (*it)->myPoint[0] * normal[0] + (*it)->myPoint[1] * normal[1] + (*it)->myPoint[2] * normal[2];
		if (valueToCheckForPlane < d && Z3i::l2Metric((*it)->myPoint, center) <= distanceMax)
			inferiorPoints.insert((*it));
	}
	return inferiorPoints;
}

template <typename WeightedPoint>
void markDifferenceBetweenPlanes(set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume,
											const Z3i::RealPoint& normalPrevious, const Z3i::Point& centerPrevious,
								 const Z3i::RealPoint& normalNext, const Z3i::Point& centerNext,
								 double distanceMax = numeric_limits<double>::max()) {
	Z3i::RealPoint next = normalNext;
	if (normalPrevious.dot(normalNext) < 0) {
		next  = -normalNext;
 	}
	set<WeightedPoint*> inferiorPointsPrevious = computePointsBelowPlane(volume, normalPrevious, centerPrevious, distanceMax);
	set<WeightedPoint*> inferiorPointsNext = computePointsBelowPlane(volume, next, centerNext, distanceMax);
		
	set<WeightedPoint*> difference;
	if (inferiorPointsNext.size() > inferiorPointsPrevious.size()) {
		set_difference(inferiorPointsNext.begin(), inferiorPointsNext.end(),
					   inferiorPointsPrevious.begin(), inferiorPointsPrevious.end(),
					   inserter(difference, difference.end()));
		
	}
	else {
		set_difference(inferiorPointsPrevious.begin(), inferiorPointsPrevious.end(),
					   inferiorPointsNext.begin(), inferiorPointsNext.end(),
					   inserter(difference, difference.end()));
	}
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		if (difference.find(*it) != difference.end()) {
			(*it)->myProcessed = true;
		}
	}
}



Z3i::DigitalSet detectBranchingPointsInNeighborhood(const Z3i::DigitalSet& branchingPoints, const Z3i::DigitalSet& setVolume, 
													const Z3i::Point& current, double radius) {

	Z3i::DigitalSet aSet(branchingPoints.domain());
	Ball<Z3i::Point> ball(current, radius);

	Z3i::Point branchingPoint;
	double distanceMin = numeric_limits<double>::max();
	bool toMark = false;
	for (auto it = branchingPoints.begin(), ite = branchingPoints.end(); it != ite; ++it) {
		if (setVolume.find(*it) != setVolume.end() && ball.contains(*it) &&
			Z3i::l2Metric(*it, current) < distanceMin) {
			toMark = true;
			branchingPoint = *it;
		    distanceMin = Z3i::l2Metric(*it, current);
		}
	}

	if (toMark) {
		double ballRadius = Z3i::l2Metric(branchingPoint, current) + 1;
		ballRadius = radius;
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

template <typename VoronoiMap>
Z3i::DigitalSet extractVoronoiCell(const Z3i::DigitalSet& backgroundSet, const VoronoiMap& voronoiMap, const Z3i::Point& point) {
	Z3i::DigitalSet aSet(backgroundSet.domain());
	for (auto it = backgroundSet.begin(), ite = backgroundSet.end(); it != ite; ++it) {
		if (voronoiMap(*it) == point) {
			aSet.insert(*it);
		}
	}
	return aSet;
}


double averageDistanceToClosestPoint(const Z3i::DigitalSet& aSet, const Z3i::DigitalSet& comparisonSet) {
	double sum = 0;
	for (auto it = aSet.begin(), ite = aSet.end(); it != ite; ++it) {
		Z3i::Point closestPoint = *min_element(comparisonSet.begin(), comparisonSet.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
				return Z3i::l2Metric(one, *it) < Z3i::l2Metric(two, *it);
			});
		sum += Z3i::l2Metric(*it, closestPoint);
	}
	return (sum / aSet.size());
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

Z3i::DigitalSet linkEachPair(const Z3i::DigitalSet& toProcess, const Z3i::Object26_6& branches, const Z3i::DigitalSet& branchingPoints) {
	Z3i::DigitalSet junctions(branches.domain());
	for (auto it = toProcess.begin(), ite = toProcess.end(); it != ite; ++it) {
		double distanceMax = 0;
		vector<Z3i::Point> resultingLine;
		for (auto itQ = toProcess.begin(), itQe = toProcess.end(); itQ != itQe; ++itQ) {
			if (*it != *itQ && !(PointUtil::areNeighbors(*it, *itQ))) {
				vector<Z3i::Point> line = PointUtil::linkTwoPoints(*it, *itQ);
				bool add = true;
				for (auto it = line.begin(), ite = line.end(); it != ite; ++it) {
					if (branches.pointSet().find(*it) == branches.end())
						add = false;
				}
				Z3i::DigitalSet lineSet(branches.domain());
				lineSet.insert(line.begin(), line.end());
			    double averageDistanceToBranch = averageDistanceToClosestPoint(lineSet, branchingPoints);
				if (averageDistanceToBranch > distanceMax && add) {
					distanceMax = averageDistanceToBranch;
					resultingLine = line;
				}
			}
		}
		junctions.insert(resultingLine.begin(), resultingLine.end());		
	}
	return junctions;
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

template <typename VCM, typename KernelFunction, typename VCMSurface>
Z3i::DigitalSet connectDisconnectedComponents(const Z3i::Object26_6& volume, const Z3i::DigitalSet& skeleton, const Z3i::DigitalSet& junctions, const VCMSurface& vcmSurface, const Z3i::DigitalSet& branchingPoints) {
	Z3i::DigitalSet lineSet(skeleton.domain());
	Z3i::Object26_6 obj(Z3i::dt26_6, junctions);
	vector<Z3i::Object26_6> objects;
	back_insert_iterator< std::vector<Z3i::Object26_6> > inserter( objects );
	unsigned int nbConnectedComponents = obj.writeComponents(inserter);
	int i = 0;
	for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it, ++i) {
		trace.progressBar(i, objects.size());
		Z3i::DigitalSet toProcess = neighborsOfSet(it->pointSet(), skeleton);
		Z3i::DigitalSet line = linkEachPair(toProcess, volume, branchingPoints);
		if (line.size() > 0) {
			Z3i::Point nearestBranchP = nearestBranchPoint(line, branchingPoints);
			Z3i::Point pointInLine = *(max_element(line.begin(), line.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
						return computeCurvature<VCM, KernelFunction>(line, one) < computeCurvature<VCM, KernelFunction>(line, two);
					}));
			double curvature = computeCurvature<VCM, KernelFunction>(line, pointInLine);
			//Change with vcm surface
			auto lambda = (vcmSurface.mapPoint2ChiVCM()).at(nearestBranchP).values;
			double curvatureAtPoint = VCMUtil::computeCurvatureJunction(lambda);
			trace.info() << curvature << " " << curvatureAtPoint << endl;
			lineSet.insert(line.begin(), line.end());
		}
	}
	return lineSet;
}

///////////////////////////////////////////////////////////////////////////////
int main( int  argc, char**  argv )
{
	QApplication application(argc,argv);

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
	typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;
	typedef WeightedPoint<Z3i::RealPoint> WeightedRealPoint;
	typedef WeightedPoint<Z3i::Point> WeightedPoint;
	typedef MetricAdjacency<Space, 3> MetricAdjacency;
	typedef WeightedPointCount<Point> WeightedPointCount;
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef DGtal::DistanceTransformation<Space, ThresholdedImage, Z3i::L2Metric> DTL2;
	typedef functors::ReverseHatPointFunction<Point, double> KernelFunction;

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
		("delta,d", po::value<int>()->default_value(1), "delta for ball radius")
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
	int delta = vm["delta"].as<int>();
	double thresholdFeature = vm["thresholdFeature"].as<double>();
	bool isDT = vm["skeleton"].as<bool>();

	Viewer3D<> viewer;
	viewer.show();
	
	Image volume = VolReader<Image>::importVol(inputFilename);
	Z3i::Domain domainVolume = volume.domain();
	Z3i::DigitalSet setVolume(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume, 
												  thresholdMin-1, thresholdMax);
//	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);
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


   
	//Construct VCM surface
	 Metric l2;
	KSpace ks;
	ks.init( volume.domain().lowerBound(),
			 volume.domain().upperBound(), true );
	SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
	Surfel bel = Surfaces<KSpace>::findABel( ks, backgroundPredicate, 1000000 );
	DigitalSurfaceContainer* container =
		new DigitalSurfaceContainer( ks, backgroundPredicate, surfAdj, bel, false  );
	DigitalSurface< DigitalSurfaceContainer > surface( container ); //acquired

	//! [DVCM3D-instantiation]
	Surfel2PointEmbedding embType = Pointels; // Could be Pointels|InnerSpel|OuterSpel;
	KernelFunction chiSurface( 1.0, r );             // hat function with support of radius r
	VCMOnSurface* vcm_surface = new VCMOnSurface( surface, embType, R, r,
					 chiSurface, dt, delta, l2, true);

	
	vector<double> curvaturePoints = extractCurvatureOnPoints(*vcm_surface);
	double threshold = otsuThreshold(curvaturePoints);
	Z3i::DigitalSet branchingPoints = computeBranchingPartsWithVCMFeature(*vcm_surface, domainVolume, threshold);
	trace.info() << threshold << endl;
	NotPointPredicate notBranching(branchingPoints);
	Z3i::Object26_6 obj(Z3i::dt26_6, branchingPoints);
	vector<Z3i::Object26_6> objects;
	back_insert_iterator< std::vector<Z3i::Object26_6> > inserter( objects );
	unsigned int nbConnectedComponents = obj.writeComponents(inserter);
	Z3i::DigitalSet maxCurvaturePoints(domainVolume);
	Matrix vcmrB, evecB;
	Z3i::RealVector evalB;
	for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
		double ratioMax = 0;
		Point maximizingCurvaturePoint;
		for (auto itPoint = it->pointSet().begin(), itPointE = it->pointSet().end(); itPoint != itPointE; ++itPoint) {
			auto lambda = (vcm_surface->mapPoint2ChiVCM()).at(*itPoint).values;
			double ratio = VCMUtil::computeCurvatureJunction(lambda); 
			if (ratio > ratioMax) {
				ratioMax = ratio;
				maximizingCurvaturePoint = *itPoint;
			}
		}
   
 		maxCurvaturePoints.insert(maximizingCurvaturePoint);
	}

	
	// map<Point, double> junctionPointToRadius;
	// for (auto it = maxCurvaturePoints.begin(), ite = maxCurvaturePoints.end(); it != ite; ++it) {
	// 	junctionPointToRadius[*it] = numeric_limits<double>::max();
	// }
	
//	delete vcm_surface;
	
	 //  for (auto it = maxCurvaturePoints.begin(), ite = maxCurvaturePoints.end(); it != ite; ++it) {
	 // 	viewer << CustomColors3D(Color::Blue, Color::Blue) << *it;
	 // }
	 // for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
	 //   	viewer << CustomColors3D(Color(0,0,50,50), Color(0,0,50,50)) << *it;
	 // }
	 // viewer << Viewer3D<>::updateDisplay;
	 // application.exec();
	 // return 0;
	
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
	Z3i::RealPoint previousNormal=normal;

	Z3i::DigitalSet connectedComponent3D(domainVolume);
	Z3i::DigitalSet branchingParts(domainVolume);
	Z3i::RealPoint realCenter;
	Z3i::Point centerOfMass;
	Z3i::Point previousCenter;
	Z3i::DigitalSet branches(setVolume.domain());

	
	trace.beginBlock("Computing skeleton");
	//Main loop to compute skeleton (stop when no vol points left to process)
	while (numberLeft > 0)
	{
		trace.progressBar((setVolumeWeighted.size() - numberLeft), setVolumeWeighted.size());		
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
				chi = KernelFunction( 1.0, radius);
			}
		}
		
		// Compute discrete plane
		connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
															 currentPoint->myPoint, normal,
															 0, radius, distanceMax, true);

	    realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);
		
		double radiusIntersection = computeRadiusFromIntersection3D<MatrixXd>(connectedComponent3D);		
		double radiusJunction = radiusForJunction(branchingPoints, *vcm_surface, realCenter, radiusIntersection);
		//Center of mass computation
		if (realCenter != Z3i::RealPoint()) {
			centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
			Z3i::DigitalSet shell = computeShell(currentPoint->myPoint, setVolume, radiusIntersection*3, radiusIntersection*5);
			int degree = computeDegree(shell);
			//int degree = 2;
			int label = (*find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
						return (wpc->myPoint == centerOfMass);
					}))->myCount;
			bool processed = false;
			for (auto it = connectedComponent3D.begin(), ite = connectedComponent3D.end(); it != ite; ++it) {
				for (auto itS = skeletonPoints.begin(), itSe = skeletonPoints.end(); itS != itSe; ++itS)  {
					if (*itS == *it)
						processed = true;
				}
			}
			
			VCMUtil::markConnectedComponent3D(setVolumeWeighted, connectedComponent3D, 0);
			if (label != 1 && !processed && degree > 1 && Z3i::l2Metric(currentPoint->myPoint, centerOfMass) <= sqrt(3)) {			   
				//detectBranchingPointsInNeighborhood(junctionPointToRadius, maxCurvaturePoints, *vcm_surface, 
				//								    normal, centerOfMass, radius);
				Z3i::DigitalSet branch = detectBranchingPointsInNeighborhood(branchingPoints, setVolume, realCenter, radiusJunction);
				branches.insert(branch.begin(), branch.end());
				
				// VCMUtil::markConnectedComponent3D(setVolumeWeighted, branch, 1);
				
				// Branching detection	
				previousNormal = normal;
				previousCenter = centerOfMass;
				skeletonPoints.insert(centerOfMass);
				viewer << CustomColors3D(Color::Red, Color::Red) << centerOfMass;
				viewer << Viewer3D<>::updateDisplay;
			    qApp->processEvents();
			}
		}
		
		//Go to next point according to normal OR to max value in DT
		VCMUtil::trackNextPoint(currentPoint, dt, setVolumeWeighted, connectedComponent3D, centerOfMass, normal);
		numberLeft = count_if(setVolumeWeighted.begin(), setVolumeWeighted.end(),
							  [&](WeightedPointCount* wpc) {
								  return (!wpc->myProcessed);
							  });
  		i++;
	}
	trace.endBlock();
	delete vcm_surface;
	
	//Discarding points being in branching parts
	for (auto it = branches.begin(), ite = branches.end(); it != ite; ++it) {
	 	auto itToErase = skeletonPoints.find(*it);
	 	if (itToErase != skeletonPoints.end())
	 		skeletonPoints.erase(itToErase);
	} 

	for (auto it=skeletonPoints.begin(), ite = skeletonPoints.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color::Red, Color::Red) << *it;
	}

	// for (auto it = junctionPointToRadius.begin(), ite = junctionPointToRadius.end(); it != ite; ++it) {
	// 	if (it->second != numeric_limits<double>::max()) {
	// 		double radiusVCM = (vcm_surface->mapPoint2Radius()).at(it->first);
	// 		Ball<Point> ballVCM(it->first, radiusVCM);
	// 		vector<Point> intersectionBallVCM = ballVCM.surfaceIntersection(setSurface);
	// 		Z3i::RealPoint normalPlaneIntersection = Statistics::computeNormalFromLinearRegression<Z3i::RealPoint>(intersectionBallVCM);
	// 		if (intersectionBallVCM.size() > 0) {
	// 			Point pointInPlane = *(intersectionBallVCM.begin());
				
	// 			// if (VCMUtil::abovePlane(it->first, normalPlaneIntersection, pointInPlane)) {
	// 			// 	normalPlaneIntersection = -normalPlaneIntersection;
	// 			// }
	// 			viewer.addLine(pointInPlane, pointInPlane+normalPlaneIntersection*10);
	// 			Ball<Point> ball(it->first, it->second);
	// 			vector<Point> pointsInBall = ball.pointsInBall();
	// 			vector<Point> restrictedPointsInBall;
	// 			for (auto it = pointsInBall.begin(), ite = pointsInBall.end(); it != ite; ++it) {
	// 				if (VCMUtil::abovePlane(*it, normalPlaneIntersection, pointInPlane)) {
	// 					restrictedPointsInBall.push_back(*it);
	// 				}
	// 			}
	// 			for (auto itB = restrictedPointsInBall.begin(), itBe = restrictedPointsInBall.end(); itB != itBe; ++itB)  {
	// 				if (setVolume.find(*itB) != setVolume.end())
	// 					branches.insert(*itB);
	// 			}
	// 		}

	// 	}
	// }
	
	for (auto it = branches.begin(), ite = branches.end(); it != ite; ++it) {
	 	viewer << CustomColors3D(Color(0,50,0,50), Color(0,50,0,50)) <<*it;
	 }
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
	DGtal::imageFromRangeAndValue(skeletonPoints.begin(), skeletonPoints.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
