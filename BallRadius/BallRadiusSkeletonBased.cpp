#include <iostream>
#include <fstream>
#include <cstdio>
#include <sstream>
#include <set>
#include <QtGui/qapplication.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/images/ImageSelector.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/geometry/surfaces/ChordGenericStandardPlaneComputer.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/shapes/parametric/Ball3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/geometry/curves/GridCurve.h"

#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"


#include "SurfaceUtils.h"
#include "SurfacePoint.h"
#include "../CreateCurve/Distance.h"
#include "../CreateCurve/Ball.h"
#include "PointUtil.h"
#include "../Tangent3D/SliceUtils.h"

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

template <typename Vector, typename Point, typename Visitor, typename Vertex>
vector<Point> createShortestPath(Visitor& visitor, const Vertex& bel, const Point& center, const map<Vertex, Point>& surfelToPoint) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	map<Vertex, Vertex> aMapPrevious;
	vector<Point> thePath;
	
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
				
			    Vector normal = SliceUtils::computeNormalFromLinearRegression<Vector>(correspondingPoints);
				Vector normal2 = SliceUtils::computeNormalFromCovarianceMatrix<Vector>(correspondingPoints);
				set<Point> correspondingSet(correspondingPoints.begin(), correspondingPoints.end());
				bool isClosed = closedCurve(correspondingSet);
				//checks if the center point belongs to the plane
				double eq = normal[0] * center[0] +
					normal[1] * center[1] +
					normal[2] * center[2];
				if (isClosed && normal != Vector::zero && eq >= -0.05 && eq <= 0.05) {
					trace.info() << normal << endl;
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
	typedef  DistanceTransformation<Z3i::Space, EqualPredicate, Z3i::L2Metric> DTL2;
	typedef Z3i::KSpace KSpace;
	typedef LightImplicitDigitalSurface<KSpace, Z3i::DigitalSet >   MyDigitalSurfaceContainer;
	typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
	typedef SurfacePoint<Z3i::Point, SCell> SurfacePoint;
	
	typedef BreadthFirstVisitor<MyDigitalSurface> MyBreadthFirstVisitor;
	typedef MyBreadthFirstVisitor::Node MyNode;
	typedef MyBreadthFirstVisitor::Size MySize;
	
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

	QApplication app(argc, argv);
	Viewer3D<> viewer;
	Cell dummy;
	viewer << SetMode3D( dummy.className(), "Basic" );
	viewer.show();

	Image img  = VolReader<Image>::importVol( inputFilename );
	Image skeleton = VolReader<Image>::importVol( inputSkeletonName );
	Z3i::Domain domain = img.domain();
	KSpace ks;
	ks.init( domain.lowerBound(), 
			 domain.upperBound(), true );
	Z3i::DigitalSet set3d (domain);
	SetFromImage<Z3i::DigitalSet>::append<Image> (set3d, img, 
												  thresholdMin, thresholdMax);
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 100000 );
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

	Color color(100,100,140,128);
	int  i = 0;
	trace.beginBlock("Computing distance map");
	for (auto it = skeleton.domain().begin(), ite = skeleton.domain().end(); it != ite; ++it) {
		if (skeleton(*it) > thresholdMin) {
			Point point = pointCorrespondingToMinDistance<Point>(*it, surfacePointSet);
			SCell surfel = convertToSurfel(point, surfaceVoxelSet);
			if (surfel != SCell()) {
				MyBreadthFirstVisitor visitor( digSurf, surfel );
				vector<Point> path = createShortestPath<Vector>(visitor, surfel, *it, surfaceVoxelSet);
				trace.info() << "found one " << endl;
				for (auto it = path.begin(), ite = path.end(); it != ite; ++it) {
					viewer << CustomColors3D(Color::Red, Color::Red) << *it;
				}
				break;
			}
			
		}
		
	
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
