#include <iostream>
#include <cstdio>
#include <set>
#include <QtGui/qapplication.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/images/ImageSelector.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/shapes/parametric/Ball3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"


#include "SurfaceUtils.h"
#include "SurfacePoint.h"
#include "../CreateCurve/Distance.h"
#include "../CreateCurve/Ball.h"
#include "PointUtil.h"
#include "../Tangent3D/SliceUtils.h"
#include "Path.h"

using namespace std;
using namespace DGtal;
using namespace Z3i;
namespace po = boost::program_options;

/**
 * Returns the point contained in surfacePointSet closest to pSkeleton
 */
template <typename Point>
Point pointCorrespondingToMinDistance(const Point& pSkeleton, const set<Point>& surfacePointSet) {
	Point minimumDistPoint = *(std::min_element(surfacePointSet.begin(), surfacePointSet.end(), [&](const Point& one, const Point& two) {
				return euclideanDistance(one, pSkeleton) < euclideanDistance(two, pSkeleton);
			}));
	return minimumDistPoint;
}

/**
 * Checks if curve is closed (all points have 2 neighbors)
 */
template <typename Point>
bool closedCurve(const set<Point>& points) {
	int cpt = 0;
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		int previouscpt = cpt;
		for (auto secondIt = points.begin(), secondIte = points.end(); secondIt != secondIte; ++secondIt) {
			if (PointUtil::areAlmostSimilar(*it, *secondIt) && *secondIt != *it) {
				cpt++;
			}
		}
		if (cpt - previouscpt < 2 || cpt - previouscpt > 3) cpt = 0;
	}
	return cpt >= points.size() * 2;
}

template <typename Point>
set<Point> createPathClosestToBall(const Ball<Point>& ball, const set<Point>& surfacePointSet, const Point & point) {
	double distance = 0;
	set<Point> curve;
	curve.insert(point);
	bool closed = false;
	while (!closed) {
		for (auto itSB = surfacePointSet.begin(), itSE = surfacePointSet.end(); itSB != itSE; ++itSB) {
			auto ballPoints = ball.pointsSurfaceBall();
			for (auto itB = ballPoints.begin(), itE = ballPoints.end(); itB != itE; ++itB) {
				if (euclideanDistance(*itB, *itSB) < distance) {
					for (const Point& pCurve : curve) {
						if (*itSB != pCurve && PointUtil::areAlmostSimilar(pCurve, *itSB))	{
							curve.insert(*itSB);
							break;
						}
					}
				}
			}
		}
		distance++;
		closed = closedCurve(curve);
	}
	return curve;
}

template <typename Path, typename Visitor, typename Vertex>
vector<Path> computeSystemOfLoops(Visitor& visitor, const Vertex& bel) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	map<Vertex, Vertex> aMapPrevious;
	vector<Path> systemOfLoops;
	
	while (!visitor.finished()) {
		node = visitor.current();
		vector<Vertex> neighbors;
		back_insert_iterator<vector<Vertex>> iter(neighbors);
		visitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			auto itMapExisting = aMapPrevious.find(*it);
			if (itMapExisting == aMapPrevious.end()) {
				aMapPrevious[*it] = node.first;
			}
			else {
				Vertex tmp = node.first;
			    Vertex tmp2 = aMapPrevious[*it];
				set<Vertex> path1, path2;			   
				while (tmp != bel) {
					path1.insert(tmp);
					tmp = aMapPrevious[tmp];
				}
				while (tmp2 != bel) {
					path2.insert(tmp2);
					tmp2 = aMapPrevious[tmp2];
				}
				set<Vertex> intersect;
				set_intersection(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(intersect, intersect.begin()));
				if (intersect.size() != 0) continue;
				vector<Vertex> correspondingPath;
				std::copy (path1.begin(), path1.end(), std::back_inserter(correspondingPath));
				std::copy (path2.begin(), path2.end(), std::back_inserter(correspondingPath));
				Path path(correspondingPath, *it);
				systemOfLoops.push_back(path);
			}
		}
		if (!visitor.finished())
			visitor.expand();
	}
	return systemOfLoops;
}

template <typename Point>
map<Point, Point> computeDSSOnPath(const vector<Point>& path) {
	typedef StandardDSS6Computer<typename vector<Point>::const_iterator,int,8> SegmentComputer;
	typedef GreedySegmentation<SegmentComputer> Segmentation;

	map<Point, Point> mapPointToDirection;
	SegmentComputer algo;
	Segmentation s(path.begin(), path.end(), algo);
	
	for (auto it = s.begin(), itE = s.end(); it!=itE; ++it) {
		SegmentComputer current(*it);
	    Point direction;
		RealPoint intercept, thickness;
		current.getParameters(direction, intercept, thickness);
		if (direction != Point::zero) {
			for (const Point& pointInDSS : current) {
				mapPointToDirection[pointInDSS] = direction;
			}
		}
	}
	return mapPointToDirection;
}

template <typename Point>
bool sameSign(const vector<Point>& points) {
	int negcptX = 0, negcptY = 0, negcptZ = 0;
	int poscptX = 0, poscptY = 0, poscptZ = 0;
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		negcptX = ((*it)[0] < 0 ) ? negcptX + 1 : negcptX;
		negcptY = ((*it)[1] < 0 ) ? negcptY + 1 : negcptY;
		negcptZ = ((*it)[2] < 0 ) ? negcptZ + 1 : negcptZ;

		poscptX = ((*it)[0] > 0 ) ? poscptX + 1 : poscptX;
		poscptY = ((*it)[1] > 0 ) ? poscptY + 1 : poscptY;
		poscptZ = ((*it)[2] > 0 ) ? poscptZ + 1 : poscptZ;
	}
	bool x = (negcptX == points.size() || poscptX == points.size());
	bool y = (negcptY == points.size() || poscptY == points.size());
	bool z = (negcptZ == points.size() || poscptZ == points.size());
	return ( (x && y) || (x && z) || (y && z));
}

template <typename Point>
bool consistentCrossProductsAlong(const vector<Point>& path, const Point& center) {
	map<Point, Point> pointToDirection = computeDSSOnPath(path);
	vector<Point> crossProducts;
	for (auto it = pointToDirection.begin(), ite = pointToDirection.end(); it != ite; ++it) {
		Point vectorToCenter = center - it->first;
		Point vectorToNext = it->second;
		Point crossProduct = vectorToNext.crossProduct(vectorToCenter);
		crossProducts.push_back(crossProduct);
	}
	return sameSign(crossProducts);
}

template <typename Surfel, typename Point>
vector<Point> selectGeodesicLoops(const vector<Path<Surfel>>& systemOfLoops, const map<Surfel, Point>& surfelToPoint, const Point& center) {
	vector<Point> geodesicLoops;
	for (auto it = systemOfLoops.begin(), ite = systemOfLoops.end(); it != ite; ++it) {
		trace.beginBlock("New path");
		vector<Point> pathPoints;
		for (auto itPath = it->begin(), itPathE = it->end(); itPath != itPathE; ++itPath) {
			pathPoints.push_back(surfelToPoint.at(*itPath));
		}
		bool consistentCrossP = consistentCrossProductsAlong(pathPoints, center);
		if (consistentCrossP)
			std::copy (pathPoints.begin(), pathPoints.end(), std::back_inserter(geodesicLoops));
		trace.endBlock();
	}
	return geodesicLoops;
}


template <typename Vector, typename Point, typename Visitor, typename Vertex>
vector<Point> createShortestPath(Visitor& visitor, const Vertex& bel, const Point& center, const map<Vertex, Point>& surfelToPoint) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	map<Vertex, Vertex> aMapPrevious;
	vector<Point> thePath;
	
	while (!visitor.finished()) {
		if (node.second == 100) break;
		node = visitor.current();
		vector<Vertex> neighbors;
		back_insert_iterator<vector<Vertex>> iter(neighbors);
		visitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			auto itMapExisting = aMapPrevious.find(*it);
			if (itMapExisting == aMapPrevious.end()) {
				aMapPrevious[*it] = node.first;
			}
			else {
			    Vertex tmp = node.first;
			    Vertex tmp2 = aMapPrevious[*it];
				vector<Point> correspondingPoints;
				set<Vertex> path1, path2;
				correspondingPoints.push_back(surfelToPoint.at(*it));
				while (tmp != bel) {
					path1.insert(tmp);
					correspondingPoints.push_back(surfelToPoint.at(tmp));
					tmp = aMapPrevious[tmp];
				}
				while (tmp2 != bel) {
					path2.insert(tmp2);
					correspondingPoints.push_back(surfelToPoint.at(tmp2));
					tmp2 = aMapPrevious[tmp2];
				}
				set<Vertex> intersect;
				set_intersection(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(intersect, intersect.begin()));
				if (intersect.size() != 0) continue;
			    
				//set<Point> correspondingSet(correspondingPoints.begin(), correspondingPoints.end());
				//bool isClosed = closedCurve(correspondingSet);
				//checks if the center point belongs to the plane
				//double eq = normal[0] * center[0] +
				//	normal[1] * center[1] +
				//	normal[2] * center[2];
				Point correspondingBel = surfelToPoint.at(bel);
				Point correspondingIntersection = surfelToPoint.at(*it);
				Point vectorToCenter = center - correspondingBel;
				Point vectorToIntersection = center + vectorToCenter;
				Point middlePoint = (correspondingIntersection + correspondingBel) / 2;
				if (/*isClosed  && normal != Vector::zero &&*/PointUtil::areAlmostSimilar(correspondingIntersection, vectorToIntersection)) {
					visitor.terminate();		  
					thePath = correspondingPoints;
				}
			}
		}
		if (!visitor.finished())
			visitor.expand();
	}
	return thePath;
}

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}


/**
 * Returns a surfel associated with a Z3 point, using a container mapping surfels to points
 */
template <typename Surfel, typename Point>
Surfel convertToSurfel(const Point& point, const map<Surfel, Point> & surfelToPoint) {
	Surfel surfel;
	for (auto it = surfelToPoint.begin(), ite = surfelToPoint.end(); it != ite; ++it) {
		if (it->second == point) {
			return it->first;
		}
	}
	return surfel;
}

int main( int argc, char** argv )
{
	typedef PointVector<3, double> Vector;
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef ImageSelector<Z3i::Domain, float>::Type FloatImage;
	typedef functors::IntervalForegroundPredicate<Image> Binarizer;
	typedef functors::EqualPointPredicate<Point> EqualPredicate;
	typedef functors::NotPointPredicate<Binarizer> NotPredicate;
	typedef Z3i::KSpace KSpace;
	typedef LightImplicitDigitalSurface<KSpace, Z3i::DigitalSet >   MyDigitalSurfaceContainer;
	typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
	typedef SurfacePoint<Z3i::Point, SCell> SurfacePoint;
	
	typedef BreadthFirstVisitor<MyDigitalSurface> MyBreadthFirstVisitor;
	typedef MyBreadthFirstVisitor::Node MyNode;
	typedef MyBreadthFirstVisitor::Size MySize;
	typedef Path<SCell> Path;
	
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file")
		("skeleton,s", po::value<std::string>(), "vol file (skeleton)")
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
	string inputFilename = vm["input"].as<std::string>();
	string inputSkeletonName = vm["skeleton"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();



	Image img  = VolReader<Image>::importVol( inputFilename );
	Image skeleton = VolReader<Image>::importVol( inputSkeletonName );
	Z3i::Domain domain = img.domain();
	KSpace ks;
	ks.init( domain.lowerBound(), 
			 domain.upperBound(), true );

	QApplication app(argc, argv);
	Viewer3D<> viewer(ks);

	Cell dummy;
	viewer << SetMode3D( dummy.className(), "Basic" );
	viewer.show();
	
	Z3i::DigitalSet set3d (domain);
	SetFromImage<Z3i::DigitalSet>::append<Image> (set3d, img, 
												  thresholdMin, thresholdMax);
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 10000000 );
	typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
	MySurfelAdjacency surfAdj( true );
	MyDigitalSurfaceContainer* ptrSurfContainer = 
		new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
	MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
	set<SCell> surfelSet;
	map<SCell, Point> surfaceVoxelSet;
	set<Point> surfacePointSet;
	for (auto it = digSurf.begin(), itE = digSurf.end(); it!=itE; ++it) {
		surfelSet.insert(*it);
	}
	
	SurfaceUtils::computeSurfelWeight<SurfacePoint>(ks, set3d, surfelSet, surfacePointSet, surfaceVoxelSet);


	Z3i::DigitalSet setskel (skeleton.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image> (setskel, skeleton, 
												  thresholdMin, thresholdMax);
	Color color(100,100,140,128);
    int number = 1;
	int i = 0;
	trace.beginBlock("Computing distance map");
	while (i < number) {
		Domain domainskel = skeleton.domain();
		auto it = select_randomly(setskel.begin(), setskel.end());	 
		Point point = pointCorrespondingToMinDistance<Point>(*it, surfacePointSet);
		SCell surfel = convertToSurfel(point, surfaceVoxelSet);
		viewer << CustomColors3D(Color::Green, Color::Green) << surfel;
		if (surfel != SCell()) {
			MyBreadthFirstVisitor visitor( digSurf, surfel );
			vector<Path> systemOfLoops = computeSystemOfLoops<Path>(visitor, surfel);
			vector<Point> geodesicLoops = selectGeodesicLoops(systemOfLoops, surfaceVoxelSet, *it);
			for (auto it = geodesicLoops.begin(), ite = geodesicLoops.end(); it != ite; ++it) {
				viewer << CustomColors3D(Color::Red, Color::Red) << *it;
			}
//			vector<Point> path = createShortestPath<Vector>(visitor, surfel, *it, surfaceVoxelSet);
			/*for (auto it = systemOfLoops.begin(), ite = systemOfLoops.end(); it != ite; ++it) {
				for (auto itpath = it->begin(), itpathe = it->end(); itpath != itpathe; ++itpath) {
					viewer << CustomColors3D(Color::Red, Color::Red) << *itpath;
				}
				}*/
		}
		i++;
	}
	trace.endBlock();

	trace.beginBlock("Displaying");
	for (auto it = set3d.begin(), ite = set3d.end(); it != ite; ++it) {
		viewer << CustomColors3D(color, color) << *it;
	}
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
	trace.endBlock();
	
	return 0;
}
