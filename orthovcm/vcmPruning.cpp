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
			int previous = 0;
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

template <typename DTL2>
Z3i::DigitalSet extractSurfacePointsWithDT( const Z3i::DigitalSet& setSurface,
											const DTL2& dt,
											const Z3i::Point& center, double radius) {
	Z3i::DigitalSet surfelSet(setSurface.domain());
	Z3i::Point closestPoint = *min_element(setSurface.begin(), setSurface.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
			return (Z3i::l2Metric(one, center) < Z3i::l2Metric(two, center));
		});
	Ball<Z3i::Point> ball(closestPoint, radius);
	vector<Z3i::Point> surfels = ball.intersection(setSurface);
	surfelSet.insert(surfels.begin(), surfels.end());
	return surfelSet;
}

template <typename DTL2>
Z3i::DigitalSet restrictEdge(const Z3i::DigitalSet& edge,
							 const Z3i::DigitalSet& branchingPoints,
							 const DTL2& dt) {
	std::vector<Z3i::Point> endPoints = CurveAnalyzer::findEndPoints(edge);
	Z3i::Point extremity;
	for (const Z3i::Point& e : endPoints) {
		if (branchingPoints.find(e) != branchingPoints.end())
			extremity = e;
	}
	double ballRadius = (dt.domain().isInside(extremity)) ? dt(extremity) : 2;
	ballRadius = (ballRadius < 2) ? 2 : ballRadius;
	Ball<Z3i::Point> ball(extremity, ballRadius);
//	vector<Z3i::Point> pointsInBall = ball.intersection(edge);
	Z3i::DigitalSet restrictedEdge(edge.domain());
	for (const Z3i::Point& e : edge) {
		if (!ball.contains(e))
			restrictedEdge.insert(e);
	}
//	restrictedEdge.insert(pointsInBall.begin(), pointsInBall.end());
	if (restrictedEdge.size() < 2) restrictedEdge = edge;
	return restrictedEdge;
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

	double thresholdInRadians = angleThreshold * M_PI / 180;

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

	//branchingPoints = reduceClustersToCenters(branchingPoints, existingSkeleton);



	//Display points
	// viewer << CustomColors3D(Color::Magenta, Color::Magenta) << branchingPoints;
	//  for (const Z3i::Point& p : endPoints) {
	//   	viewer << CustomColors3D(Color::Green, Color::Green) << p;
	//  }
	 // for (const Z3i::Point& p : existingSkeletonOrdered)
	 // 	 viewer << CustomColors3D(Color::Red, Color::Red) << p;
//	 viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
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
	//  }
	//  viewer << Viewer3D<>::updateDisplay;
	//  application.exec();
	//  return 0;


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



	vector<Z3i::DigitalSet> planes;
	/* Working with 3 planes */
	Z3i::DigitalSet processedEdges(existingSkeleton.domain());
	Z3i::DigitalSet skeletonPoints(setVolume.domain());



	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );


	//Construct VCM surface
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
	int maxLabel = 0;
	for (GraphEdge* graphEdge : hierarchicalGraph) {
		int currentLabel = graphEdge->getLabel();
		if (currentLabel > maxLabel && currentLabel != std::numeric_limits<int>::max())
			maxLabel = currentLabel;
	}
	trace.info() << maxLabel << endl;
	Z3i::DigitalSet newSkeleton = existingSkeleton, previousSkeleton = existingSkeleton;
	for (int label = 1; label <= maxLabel; label++) {
		Z3i::DigitalSet deletedEdges(existingSkeleton.domain());
		trace.info() << hierarchicalGraph.size() << endl;
		for (GraphEdge* graphEdge : hierarchicalGraph) {
			Z3i::DigitalSet edge = graphEdge->pointSet();
			if (edge.size() == 0 //|| graphEdge->getLabel() != 1
				) continue;
			Z3i::DigitalSet restrictedEdge = edge;//restrictEdge(edge, branchingPoints, dt);

			trace.progressBar(i, edgeGraph.size());

			VCM vcm(R, ceil(r), l2, false);
			vcm.init(restrictedEdge.begin(), restrictedEdge.end());
			bool add = true;
			//Main loop to compute skeleton (stop when no vol points left to process)
			int cpt = 0;
			Z3i::DigitalSet smoothedEdge(edge.domain());
			double sumAngle = 0;
			for (const Z3i::Point& s : restrictedEdge)
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
				radius = edge.size()*0.4;

				connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
																	 s, normal,
																	 0, radius, distanceMax, 26, true);

				//Z3i::DigitalSet surfelSet = extractRelevantSurfacePointsForVolVCM(connectedComponent3D, dt, s);

				double distanceMaxEdge = dt(*max_element(edge.begin(), edge.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
							return dt(one) < dt(two);
						}))+1;
				Z3i::DigitalSet surfelSet = extractSurfacePointsWithDT(setSurface, dt, s, distanceMaxEdge);
				vector<Z3i::Point> surfels(surfelSet.begin(), surfelSet.end());

				chiSurface = KernelFunction(1.0, distanceMaxEdge);
				vcmSurface.setMySmallR(distanceMaxEdge);


				Z3i::Point farthestPoint = *max_element(surfelSet.begin(), surfelSet.end(), [&](const Z3i::Point& one,
																								const Z3i::Point& two)  {
															return Z3i::l2Metric(one, s) < Z3i::l2Metric(s, two);
														});
				double radiusSurface =dt(s) + delta;
				Z3i::RealPoint normalSurface = VCMUtil::computeNormalFromVCM(s, vcmSurface, chiSurface, 0, Z3i::RealVector(), surfels);
				// connectedComponent3D = VCMUtil::computeDiscretePlane(vcmSurface, chiSurface, domainVolume, setVolumeWeighted,
				// 													 s, normalSurface,
				// 													 0, radiusSurface, distanceMax, 26, true);
				double angle = normalSurface.cosineSimilarity(normal);
				double otherAngle = normalSurface.cosineSimilarity(-normal);
				angle = (angle < otherAngle) ? angle : otherAngle;

				bool keep = ((angle < thresholdInRadians) // &&
							 // edge.size() >= 2 * distanceMaxEdge)
							 ||
							 surfels.size() < 3 // || distanceMaxEdge <= 3
					);
				sumAngle += angle;
				cpt++;

			}
			sumAngle /= cpt;
			if (sumAngle > thresholdInRadians) {
				//	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << edge;
				Z3i::DigitalSet difference(newSkeleton.domain());
				for (const Z3i::Point& p : newSkeleton) {
					bool found = false;
					if (edge.find(p) != edge.end() &&
						branchingPoints.find(p) == branchingPoints.end())
						found = true;
					if (!found) {
						difference.insert(p);
					}
				}
				Z3i::Object26_6 objDiff(Z3i::dt26_6, difference);
				vector<Z3i::Object26_6> objects;
				back_insert_iterator<vector<Z3i::Object26_6>> inserter(objects);
				unsigned int nbCC = objDiff.writeComponents(inserter);
				if (nbCC == 1
					) {
				    newSkeleton = (objects.begin())->pointSet();

				}
			}

			i++;
		}

	    endPointsV = CurveAnalyzer::findEndPoints(newSkeleton);
		p = (*endPointsV.begin());

	    branchingPoints = Z3i::DigitalSet(domainVolume);
	    existingSkeletonOrdered = CurveAnalyzer::curveTraversalForGraphDecomposition(branchingPoints,
																					 newSkeleton,
																					 p);
	    endPoints = vector<Z3i::Point>();
		for (const Z3i::Point& p : endPointsV) {
			if (branchingPoints.find(p) == branchingPoints.end())
				endPoints.push_back(p);
		}
	    edgeGraph = CurveAnalyzer::constructGraph(existingSkeletonOrdered, branchingPoints);
		hierarchicalGraph = CurveAnalyzer::hierarchicalDecomposition(edgeGraph, endPoints, branchingPoints);

	}

	trace.endBlock();
	viewer << CustomColors3D(Color::Red, Color::Red) << newSkeleton;

	//second pass
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
	    viewer << CustomColors3D(Color(0,0,30,10), Color(0,0,30,10)) << (*it)->myPoint;
	}

	Image outImage2(volume.domain());
	DGtal::imageFromRangeAndValue(newSkeleton.begin(), newSkeleton.end(), outImage2, 10);
	VolWriter<Image>::exportVol(outFilename, outImage2);
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
