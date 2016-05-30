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



using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

class GraphEdge {
public:
	GraphEdge(const Z3i::DigitalSet& points, int label) : myPoints(points), myLabel(label) {}
	int getLabel() const { return myLabel; }
	void setLabel(int label) { myLabel = label; }
	Z3i::DigitalSet pointSet() const { return myPoints; }
private:
	Z3i::DigitalSet myPoints;
	int myLabel;
};

class Concatenation {
public:
	Concatenation(const vector<Z3i::DigitalSet>& edges, int level) : myEdges(edges), myLevel(level) {}

	double computeAverageFunction(const std::function<double(const Z3i::DigitalSet& aSet)>& func,
								  const std::function<bool(const Z3i::DigitalSet& aSet)>& pred = {}) const {
		double sumValue = 0;
		int cpt = 0;
		for (const Z3i::DigitalSet& edge : myEdges) {
			bool checkPred = true;
			if (pred) checkPred = pred(edge);
			if (checkPred) {
				sumValue += func(edge);
				cpt++;
			}
		}
	    if (cpt == 0) return 0;
		return sumValue / cpt;
	}
	
public:
	vector<Z3i::DigitalSet> myEdges;
	int myLevel;
};

class LevelConcatenation {
public:
	LevelConcatenation(const vector<Concatenation>& concats, int level) : myConcatenations(concats), myLevel(level) {}

	double computeAverageLevelFunction(const std::function<double(const Z3i::DigitalSet& aSet)>& func,
									   const std::function<bool(const Z3i::DigitalSet& aSet)>& pred = {}) const {
		double sumValue = 0;
		int cpt = 0;
		for (const Concatenation& concat : myConcatenations) {
			sumValue += concat.computeAverageFunction(func, pred);
		}
		return sumValue / myConcatenations.size();
	}

public:
	vector<Concatenation> myConcatenations;
	int myLevel;
};

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

Z3i::DigitalSet detectCriticalPoints(const Z3i::DigitalSet& skeleton) {
	Z3i::Object26_6 obj(Z3i::dt26_6, skeleton);
	Z3i::DigitalSet criticalPoints(skeleton.domain());
	for (const Z3i::Point& s : skeleton) {
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
		obj.writeNeighbors(inserter, s);
		if (neighbors.size() > 2) {
			criticalPoints.insert(s);
		}
	}
	return criticalPoints;
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

double computeAngle(const Z3i::Point& b,
			 const Z3i::Point& m,
			 const Z3i::Point& e) {
	Z3i::RealVector v1 = (b - m).getNormalized();
	Z3i::RealVector v2 = (e - m).getNormalized();
    double dot = v1.dot(v2);
	double angle = acos(dot);
	return angle;
}

template <typename Segmentation, typename Domain>
Z3i::DigitalSet dominantPointDetection(const Segmentation& segmentation,
									   const vector<Z3i::Point>& skeletonOrdered,
									   const Domain& domain) {

	typedef typename Segmentation::SegmentComputerIterator DSS;
	
	Z3i::DigitalSet dominantPoints(domain);	
	DSS q = (segmentation.begin());
	DSS p = (++segmentation.begin());
	while (p != segmentation.end()) {
		Z3i::Point endQ = *(--(q->end()));
		Z3i::Point beginP = *(p->begin());
		int posP = find(skeletonOrdered.begin(), skeletonOrdered.end(), beginP) - skeletonOrdered.begin();
		int posQ = find(skeletonOrdered.begin(), skeletonOrdered.end(), endQ) - skeletonOrdered.begin();
	    if (posP < posQ) {
			Z3i::Point beginQ = *(q->begin());
			Z3i::Point endP = *(--(p->end()));
			double min = numeric_limits<double>::max();
			Z3i::Point dominantCandidate;
			for (int i = posP; i < posQ; i++) {
				Z3i::Point pointI = skeletonOrdered[i];
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
Z3i::DigitalSet branchingPointDetection(const Segmentation& segmentation,
										const vector<Z3i::Point>& skeletonOrdered,
										const Domain& domain) {

	typedef typename Segmentation::SegmentComputerIterator DSSIterator;
	typedef typename Segmentation::SegmentComputer DSS; 

	Z3i::DigitalSet branchingPoints(domain);	

	for (const Z3i::Point& s : skeletonOrdered) {

		vector<DSS> segments;
		
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
				Z3i::Point beginI = *(si.begin());
				Z3i::Point endI = *(--(si.end()));
				Z3i::Point beginJ = *(sj.begin());
				Z3i::Point endJ = *(--(sj.end()));
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


vector<Z3i::DigitalSet> constructGraph(const vector<Z3i::Point>& orderedCurve,
									   const Z3i::DigitalSet& constraint) {
	vector<Z3i::DigitalSet> graph;
	int index = 0;
	graph.push_back(Z3i::DigitalSet(constraint.domain()));
	for (int i = 0, end = orderedCurve.size(); i < end; i++) {
		Z3i::Point current = orderedCurve[i];		
		if (constraint.find(current) != constraint.end()) {
			index++;
			graph.push_back(Z3i::DigitalSet(constraint.domain()));
			graph[index].insert(current);
		}
		graph[index].insert(current);
	}
	return graph;
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
	double delta = dt(b) - dt(end);
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

bool sameSet(const Z3i::DigitalSet& first, const Z3i::DigitalSet& second) {
	if (first.size() != second.size()) return false;
	unsigned int cpt = 0;
	for (const Z3i::Point& f : first) {
		if (second.find(f) != second.end())
			cpt++;
	}
	return (cpt == first.size());
}

vector<GraphEdge*> neighboringEdges(const vector<GraphEdge*>& edges,
									const Z3i::DigitalSet& currentEdge,
									const Z3i::DigitalSet& branchingPoints) {
	typedef MetricAdjacency<Z3i::Space, 3> MetricAdjacency;
	
	vector<GraphEdge*> neighbors;
	Z3i::Point branchPoint;
	for (const Z3i::Point& b : branchingPoints) {
		if (currentEdge.find(b) != currentEdge.end() ) {
			branchPoint = b;
		}
	}

	vector<Z3i::Point> nb;
	back_insert_iterator<vector<Z3i::Point>> inserter(nb);
	MetricAdjacency::writeNeighbors(inserter, branchPoint);
	
	for (GraphEdge* edge : edges) {
		Z3i::DigitalSet setEdge = edge->pointSet();
		if (sameSet(setEdge, currentEdge)) continue;
		for (const Z3i::Point& n : nb) {
			if (setEdge.find(n) != setEdge.end())
				neighbors.push_back(edge);
		}
	}
	
	return neighbors;
}

vector<GraphEdge*> hierarchicalDecomposition(const vector<Z3i::DigitalSet>& edges,											
											const vector<Z3i::Point>& endPoints,
											const Z3i::DigitalSet& branchingPoints) {

	queue<GraphEdge*> edgeQueue;
	vector<GraphEdge*> hierarchyGraph;
	for (const Z3i::DigitalSet& edge : edges) {
		GraphEdge* levelEdge = new GraphEdge(edge, std::numeric_limits<int>::max());
	
		for (const Z3i::Point& e : endPoints) {
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
		Z3i::DigitalSet edgeCurrentSet = edgeCurrent->pointSet();
		edgeQueue.pop();
		vector<GraphEdge*> neighborEdges = neighboringEdges(hierarchyGraph, edgeCurrentSet, branchingPoints);
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

Z3i::DigitalSet ensureConnexity(const Z3i::DigitalSet& set) {
	Z3i::DigitalSet cleanSet(set.domain());
	Z3i::Object26_6 obj(Z3i::dt26_6, set);
	Z3i::DigitalSet & S = obj.pointSet();
	cleanSet = S;
	for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
		Z3i::Object26_6 obj(Z3i::dt26_6, cleanSet);
		if (obj.isSimple(*it)) {
		    cleanSet.erase(*it);
		}
	}

	return cleanSet;
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

	SetFromImage<Z3i::DigitalSet>::append<Image>(setSkeleton, skeleton, thresholdMin-1, thresholdMax);
	
	Z3i::DigitalSet existingSkeleton = ensureConnexity(setSkeleton);
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

	Z3i::DigitalSet branchingPoints(domainVolume);
    unsigned int previous = 0;
	while ( !visitor.finished() ) 
	{
  		node = visitor.current();
		if ( existingSkeleton.find(node.first) != existingSkeleton.end() ) { //is inside domain			
		    if (node.second <= previous) {
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
	// Z3i::DigitalSet branchingPoints = detectCriticalPoints(existingSkeleton);
	// branchingPoints = reduceClustersToCenters(branchingPoints);
	
	vector<Z3i::DigitalSet> edgeGraph = constructGraph(existingSkeletonOrdered, branchingPoints);
	vector<Z3i::Point> endPoints = findEndPoints(existingSkeleton);
	Z3i::DigitalSet endPointSet(domainSkeleton);
	endPointSet.insert(endPoints.begin(), endPoints.end());
	vector<GraphEdge*> hierarchicalGraph = hierarchicalDecomposition(edgeGraph, endPoints, branchingPoints);
																	
	//Display points
	//  for (const Z3i::Point& p : endPoints) {
	//  	viewer << CustomColors3D(Color::Green, Color::Green) << p;
	//  }
	// viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
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

	int i = 0;
   

	double distanceMax = dt(*max_element(existingSkeleton.begin(), existingSkeleton.end(), [&](const Z3i::Point& one, const Z3i::Point& two) {
				return dt(one) < dt(two);
			}));
	
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
				
				vector<GraphEdge*> neighborEdges= neighboringEdges(hierarchicalGraph, currentEdge, branchingPoints);
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
	double nLocMax1 = 0;
	for (const LevelConcatenation& levelConcat : groupConcatenations) {
		nLocMax1 = levelConcat.computeAverageLevelFunction(numberLocMaxFunction, predicate);
	}
	nLocMax1 /= groupConcatenations.size();


	Z3i::DigitalSet toPrune(domainSkeleton);
	double criterionLength = (la0 + la1) / 2.;
	for (const Concatenation& concat : levelConcat1.myConcatenations) {
		double lengthH = concat.computeAverageFunction(lengthFunction);
		double nLocMaxH = concat.computeAverageFunction(numberLocMaxFunction);
		double deltaH = concat.computeAverageFunction(deltaFunction);
		if ( deltaH > -10 && nLocMaxH < 0.3 * nLocMax1 &&
			 (lengthH < criterionLength || nLocMaxH <= nLocMax0)
			) {
			for (const Z3i::DigitalSet& aSet : concat.myEdges)
				toPrune.insert(aSet.begin(), aSet.end());			
		}
	
	}


	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << toPrune;
	viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	
	trace.endBlock();
	   
	viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
	
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
