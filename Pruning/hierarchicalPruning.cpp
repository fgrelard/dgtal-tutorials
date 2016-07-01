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
#include "surface/SurfaceTraversal.h"
#include "Statistics.h"
#include "shapes/Ball.h"
///////////////////////////////////////////////////////////////////////////////
#include "ReverseHatPointFunction.h"

#include "graph/GraphEdge.h"
#include "graph/Concatenation.h"
#include "graph/LevelConcatenation.h"
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
double deltaEdge(const DTL2& dt,
			 const Z3i::DigitalSet& edge,
			 const Z3i::DigitalSet& branch,
			 const Z3i::DigitalSet& endPoints) {

	Z3i::Point b, end;
	for (const Z3i::Point& e : edge) {
		if (branch.find(e) != branch.end())
			b = e;
		if (endPoints.find(e) != endPoints.end())
			end = e;
	}
	double delta = 0;
	if (b != Z3i::Point() && end != Z3i::Point())
		delta = dt(b) - dt(end);
	return delta;
}

template <typename DTL2>
int numberOfLocalMaximaDT(const DTL2& dt,
						  const Z3i::DigitalSet& aSet) {
	typedef MetricAdjacency<Z3i::Space, 3> MAdj;

	int localMax = 0;

	for (const Z3i::Point& p : aSet) {
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
		MAdj::writeNeighbors(inserter, p);
		double currentDTValue = dt(p);
		bool isLocalMaximum = true;
		for (const Z3i::Point& n: neighbors) {
			if (dt.domain().isInside(n) && dt(n) > currentDTValue)
				isLocalMaximum = false;
		}
		if (isLocalMaximum) localMax++;
	}
	return localMax;
}

unsigned int lengthEdge(const Z3i::DigitalSet& edge) {
	return edge.size();
}

template <typename DTL2>
double amountInformationLostRatio(const DTL2& dt,
					   const Z3i::DigitalSet& edge,
					   const Z3i::DigitalSet& branch,
					   const Z3i::DigitalSet& endPoints) {
	Z3i::Point b, end;
	for (const Z3i::Point& e : edge) {
		if (branch.find(e) != branch.end())
			b = e;
		if (endPoints.find(e) != endPoints.end())
			end = e;
	}
	double length = lengthEdge(edge);
	double delta = deltaEdge(dt, edge, branch, endPoints);
	double Rh = dt(b);
	double d = Z3i::l2Metric(b, end);
	double qratio = length * (d - delta) / Rh;
	return qratio;
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
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
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

	SetFromImage<Z3i::DigitalSet>::append<Image>(setSkeleton, skeleton, thresholdMin-1, thresholdMax);

	Z3i::DigitalSet existingSkeleton = CurveAnalyzer::ensureConnexity(setSkeleton);
	typedef StandardDSS6Computer<vector<Point>::iterator,int,8> SegmentComputer;
	typedef GreedySegmentation<SegmentComputer> Segmentation;
	vector<Point> existingSkeletonOrdered;
		vector<Z3i::Point> endPointsV = CurveAnalyzer::findEndPoints(existingSkeleton);
	Z3i::Point p = (*endPointsV.begin());
	// We have to visit the direct neighbours in order to have a container with voxels
	// ordered sequentially by their connexity
	// Otherwise we have a point container with points which are not neighbours
	// and this impairs maximal segment recognition
	typedef Z3i::Object26_6 Graph;
	typedef DepthFirstVisitor<Graph, set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;
	typedef GraphVisitorRange<Visitor> VisitorRange;
	Graph graph(Z3i::dt26_6, existingSkeleton);
	Visitor visitor( graph, p );
	MyNode node;

	Z3i::DigitalSet branchingPoints(domainVolume);
    pair<Z3i::Point, double> previous;

	while ( !visitor.finished() )
	{
  		node = visitor.current();
		if (node.second != 0 && ((int)node.second - previous.second) <= 0) {
			vector<Z3i::Point> neighbors;
			back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
			MetricAdjacency::writeNeighbors(inserter, node.first);
			double minDistance = std::numeric_limits<double>::max();
			Z3i::Point cand;
			for (const Z3i::Point& n : neighbors) {
				if (find(existingSkeletonOrdered.begin(), existingSkeletonOrdered.end(), n) != existingSkeletonOrdered.end()) {
					double currentDistance = Z3i::l2Metric(n, node.first);
					if (currentDistance < minDistance) {
						minDistance = currentDistance;
						cand = n;
					}
				}
			}
			branchingPoints.insert(cand);
			existingSkeletonOrdered.push_back(cand);

		}
		previous = node;
		existingSkeletonOrdered.push_back(node.first);
		visitor.expand();
	}
	SegmentComputer algo;
	Segmentation s(existingSkeletonOrdered.begin(), existingSkeletonOrdered.end(), algo);
	s.setMode("MostCentered++");
	typename Segmentation::SegmentComputerIterator itseg = s.begin();
	typename Segmentation::SegmentComputerIterator end = s.end();
//	Z3i::DigitalSet branchingPoints = branchingPointDetection(s, existingSkeletonOrdered, domainVolume);
	// Z3i::DigitalSet branchingPoints = detectCriticalPoints(existingSkeleton);
	// branchingPoints = reduceClustersToCenters(branchingPoints);

	vector<Z3i::DigitalSet> edgeGraph = CurveAnalyzer::constructGraph(existingSkeletonOrdered, branchingPoints);
	vector<Z3i::Point> endPoints = CurveAnalyzer::findEndPoints(existingSkeleton);
	Z3i::DigitalSet endPointSet(domainSkeleton);
	endPointSet.insert(endPoints.begin(), endPoints.end());
	vector<GraphEdge*> hierarchicalGraph = CurveAnalyzer::hierarchicalDecomposition(edgeGraph, endPoints, branchingPoints);

	//Display points
	//  for (const Z3i::Point& p : endPoints) {
	//  	viewer << CustomColors3D(Color::Green, Color::Green) << p;
	//  }
	//viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	//  for (GraphEdge* edge : hierarchicalGraph) {
	// 	 int label = edge->getLabel();
	// 	 if (label > 3) {
	// 		 int r = (label * 64) % 256;
	// 		 int g = (label* 128)%256;
	// 		 int b = (label* 192)%256;
	// 		 viewer << CustomColors3D(Color(r,g,b), Color(r,g,b)) << edge->pointSet();
	// 	 }
	// 	 else {
	// 		 viewer << CustomColors3D(Color::Red, Color::Red)  << edge->pointSet();
	// 	 }
	//  }

	// viewer << Viewer3D<>::updateDisplay;
	// application.exec();

	// for (GraphEdge* edge : hierarchicalGraph) {
	// 	if (edge->getLabel() == 1) {
	// 		viewer << CustomColors3D(Color::Yellow, Color::Yellow);
	// 	}
	// 	else
	// 		viewer << CustomColors3D(Color::Red, Color::Red);
	// 	viewer << edge->pointSet();
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



	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );

	trace.beginBlock("Pruning skeleton");
	int maxLabel = 0;
	for (GraphEdge* graphEdge : hierarchicalGraph) {
		int currentLabel = graphEdge->getLabel();
		if (currentLabel > maxLabel && currentLabel != std::numeric_limits<int>::max())
			maxLabel = currentLabel;
	}
	trace.info() << maxLabel << endl;

	std::function<double(const Z3i::DigitalSet& aSet)> deltaFunction = [&](const Z3i::DigitalSet& aSet) {
		return deltaEdge(dt, aSet, branchingPoints, endPointSet);
	};
	std::function<double(const Z3i::DigitalSet& aSet)> lengthFunction = [&](const Z3i::DigitalSet& aSet) {
		return lengthEdge(aSet);
	};

	std::function<double(const Z3i::DigitalSet& aSet)> numberLocMaxFunction = [&](const Z3i::DigitalSet& aSet) {
		return numberOfLocalMaximaDT(dt, aSet);
	};
	std::function<double(const Z3i::DigitalSet& aSet)> qRatioFunction = [&](const Z3i::DigitalSet& aSet) {
		return amountInformationLostRatio(dt, aSet, branchingPoints, endPointSet);
	};


	vector<LevelConcatenation> groupConcatenations;
	for (int level = 1; level <= maxLabel; level++) {
		vector<Concatenation> concatenations;
		for (GraphEdge* graphEdge : hierarchicalGraph) {
			if (graphEdge->getLabel() == 1) { // Extremity
				vector<Z3i::DigitalSet> edges;
				Z3i::DigitalSet currentEdge = graphEdge->pointSet();
				edges.push_back(currentEdge);

				vector<GraphEdge*> neighborEdges= CurveAnalyzer::neighboringEdges(hierarchicalGraph, currentEdge, branchingPoints);
				int localMaxLabel = 0;
				for (GraphEdge* g : neighborEdges) {
					int currentLabel = g->getLabel();
					if (currentLabel == 1) continue;
					if (currentLabel > localMaxLabel && currentLabel <= level)
						localMaxLabel = currentLabel;
				}

				for (GraphEdge* g : neighborEdges) {
					if (g->getLabel() != 1 && g->getLabel() == localMaxLabel) {
						edges.push_back(g->pointSet());
					}
				}

				Concatenation concat(edges, level);
				concatenations.push_back(concat);

			}
		}
		LevelConcatenation levelConcat(concatenations, level);
		groupConcatenations.push_back(levelConcat);
	}

	LevelConcatenation levelConcat1 = *(find_if(groupConcatenations.begin(), groupConcatenations.end(), [&](const LevelConcatenation& lConcat) {
			return lConcat.myLevel == 1;
			}));

	double la1 = levelConcat1.computeAverageLevelFunction(lengthFunction);
	std::function<bool(const Z3i::DigitalSet& aSet)> predicate = [&](const Z3i::DigitalSet& aSet) {
	 	return lengthEdge(aSet) < la1;
    };

	double la0 = levelConcat1.computeAverageLevelFunction(lengthFunction, predicate);
	double nLocMax0 = levelConcat1.computeAverageLevelFunction(numberLocMaxFunction, predicate);
	double nLocMax1 = levelConcat1.computeAverageLevelFunction(numberLocMaxFunction);

	Z3i::DigitalSet toPrune(domainSkeleton);
	double criterionLength = (la0 + la1) / 2.;
	Z3i::DigitalSet resultingSkeleton(existingSkeleton.domain());
	for (const Concatenation& concat : levelConcat1.myConcatenations) {
		double lengthH = concat.computeAverageFunction(lengthFunction);
		double nLocMaxH = concat.computeAverageFunction(numberLocMaxFunction);
		double deltaH = concat.computeAverageFunction(deltaFunction);
		if (  deltaH > -10 && nLocMaxH < 0.3 * nLocMax1 &&
			 (lengthH < criterionLength || nLocMaxH <= nLocMax0)
			) {
			for (const Z3i::DigitalSet& aSet : concat.myEdges) {
				toPrune.insert(aSet.begin(), aSet.end());
			}
		}
	}

	for (const Z3i::Point& p : existingSkeleton) {
		if (toPrune.find(p) == toPrune.end() || branchingPoints.find(p) != branchingPoints.end())
			resultingSkeleton.insert(p);
	}

	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << toPrune;
	viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;

	trace.endBlock();

	//second pass
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
	    viewer << CustomColors3D(Color(0,0,120,10), Color(0,0,120,10)) << (*it)->myPoint;
	}


	Image outImage(volume.domain());
	DGtal::imageFromRangeAndValue(resultingSkeleton.begin(), resultingSkeleton.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
