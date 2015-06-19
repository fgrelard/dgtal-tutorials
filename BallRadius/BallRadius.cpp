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
//LICENSE-END
/**
 * @file testSegmentation.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 *
 * @date 2011/07/22
 *
 * This file is part of the DGtal library
 */

/**
 * Description of testSegmentation <p>
 * Aim: simple test of \ref GreedySegmentation and \ref SaturatedSegmentation
 */

#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <iterator>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/base/Exceptions.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/base/Circulator.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <Eigen/Dense>

#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/kernel/BasicPointPredicates.h"



// Local includes
#include "geometry/RosenProffittLengthEstimator.h"
#include "surface/SurfacePoint.h"
#include "surface/WeightedDigitalSurface.h"
#include "geometry/Distance.h"
#include "geometry/SliceUtils.h"
#include "geometry/Pencil.h"
#include "geometry/MSTTangent.h"
#include "surface/RadialSurface.h"
#include "Statistics.h"
#include "shapes/Ellipse.h"
#include "surface/SurfaceUtils.h"
#include "geometry/PointUtil.h"

using namespace DGtal;
using namespace std;
using namespace Z3i;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GreedySegmentation
///////////////////////////////////////////////////////////////////////////////


template <typename Surfel>
void visualizePath(const Surfel& start, const Surfel& end, map<Surfel, Surfel>& previous, Viewer3D<>& viewer, const Color& color) {
	Surfel tmp = end;
	while (tmp != start) {
		viewer << CustomColors3D(color, color) << tmp;
	    tmp = previous[tmp];
	}

}




template <typename Ellipse>
bool checkIfEllipseFits(const Ellipse& fittedEllipse, const Ellipse& intersectedEllipse) {
	int majorFitted = (int)fittedEllipse.myMajorAxis;
	int majorIntersected = (int)intersectedEllipse.myMajorAxis;

	int minorFitted = (int)fittedEllipse.myMinorAxis;
	int minorIntersected = (int)intersectedEllipse.myMinorAxis;

	bool isMajor = majorFitted != 0 && (majorFitted == majorIntersected || majorFitted == majorIntersected - 1 || majorFitted == majorIntersected + 1);
	bool isMinor = minorFitted != 0 && (minorFitted == minorIntersected || minorFitted == minorIntersected - 1 || minorFitted == minorIntersected + 1);

	return isMajor && isMinor;
}

template <typename Domain, typename Point, typename Surfel>
vector<Point> surfacePointsOnPaths(const map<Surfel, Point>& surfelToPoint, const set<Point>& surfacePointSet, const set<Surfel>& path) {
	vector<Point> pathPoints;
	for (auto it = path.begin(), ite = path.end(); it != ite; ++it) {
		pathPoints.push_back(surfelToPoint.at(*it));
	}
	Domain domain = PointUtil::computeBoundingBox<Domain>(pathPoints);
	vector<Point> surfacePoints;
	for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
	    if (surfacePointSet.find(*it) != surfacePointSet.end()) {
			surfacePoints.push_back(*it);
		}
	}
	return surfacePoints;
}

template <typename Point>
bool isAlmostSurfacePoint(const Point& point, const set<Point>& surfacePointSet) {
	for (auto it = surfacePointSet.begin(), ite = surfacePointSet.end(); it != ite; ++it) {
		if (PointUtil::areAlmostSimilar(*it, point)) return true;
	}
	return false;
}

template <typename Domain, typename Surfel, typename Point>
bool checkSymmetry(const map<Surfel, Point>& surfelToPoint, const set<Point>& surfacePointSet, const set<Surfel>& path1, const set<Surfel>& path2, const Point & center) {
	if (path2.size() == 0) return false;
	int foundOneSymmetry = 0;
	vector<Point> pathPoints = surfacePointsOnPaths<Domain>(surfelToPoint, surfacePointSet, path2);
	vector<Point> otherPathPoints = surfacePointsOnPaths<Domain>(surfelToPoint, surfacePointSet, path1);
	vector<Point> intersection;
	sort(pathPoints.begin(), pathPoints.end());
	sort(otherPathPoints.begin(), otherPathPoints.end());
	set_intersection(pathPoints.begin(), pathPoints.end(), otherPathPoints.begin(), otherPathPoints.end(), back_inserter(intersection));
	//If the two sets of points overlap, then we didnt find the cross section corresponding to a tube (ellipsoidal)
	if (intersection.size() >= 0.5  * path1.size()) return false;
	
	int cpt = 0;
	for (auto it = path1.begin(), ite = path1.end(); it != ite; ++it) {
		Point middle = (surfelToPoint.at(*(path1.begin())) + surfelToPoint.at(*it)) /2;
		if (isAlmostSurfacePoint(middle, surfacePointSet)) {
			cpt++;
		}
	}
	//If we have all the points in the path almost being surface points, then the path found does not correspond to a cross section
	if (cpt == path1.size()) return false;
	for (auto it = path1.begin(), itE = path1.end(); it != itE; ++it) {
		Point currentPoint = surfelToPoint.at(*it);
		Point vectorToCenter = center - currentPoint;
		Point symmetryCurrent = center + vectorToCenter;
		for (auto itS = pathPoints.begin(), itSE = pathPoints.end(); itS != itSE; ++itS) {
			Point putativeSymmetric = *itS;
			if (PointUtil::areAlmostSimilar(putativeSymmetric, symmetryCurrent) && !PointUtil::areAlmostSimilar(center, symmetryCurrent) && !PointUtil::areAlmostSimilar(center, currentPoint) && !PointUtil::areAlmostSimilar(putativeSymmetric, currentPoint)) {
				foundOneSymmetry++;
				break;
			}
		}
	}
	return (foundOneSymmetry == path1.size());
}



template <typename Point>
vector<Point> nearestPointsFromCenter(const vector<Point> & points, const Point & center, double radius) {
	vector<Point> nearestPointV;
	for (const Point& point : points) {
		if (euclideanDistance(point, center) <= radius) {
			nearestPointV.push_back(point);
		}
	}
	return nearestPointV;
}

template <typename Surfel, typename Point>
bool areNeighbours(const set<Surfel>& path1, const set<Surfel>& path2, const map<Surfel, Point>& surfelToPoint) {
	int cpt = 0;
	for (const auto & surfelInPOne : path1) {
		Point pp1 = surfelToPoint.at(surfelInPOne);
		for (const auto & surfelInPTwo : path2) {
			Point pp2 = surfelToPoint.at(surfelInPTwo);
			if (areAlmostSimilar(pp1, pp2))
				cpt++;
		}
	}
	return cpt>4;
}

template <typename Vertex, typename Tracker, typename DirIterator>
vector<Vertex> computeMapInOneDirection(Vertex& s, Tracker tracker, DirIterator& itDirs, const int notDirection) {
	vector<Vertex> aVectorNeighbors;
	for (; itDirs != 0; ++itDirs) {
		if ( tracker->adjacent( s, *itDirs, true )) {
			aVectorNeighbors.push_back(s);
		}
		if ( tracker->adjacent( s, *itDirs, false )) {
			if (*itDirs == notDirection)
				continue;
			aVectorNeighbors.push_back(s);
		}
	}
	return aVectorNeighbors;
}

template <typename Point>
double estimateDSSLength(const vector<Point>& points) {
	typedef StandardDSS6Computer<vector<Z3i::Point>::const_iterator,int,8> SegmentComputer;  
	typedef GreedySegmentation<SegmentComputer> Segmentation;
	SegmentComputer algo;
	Segmentation s(points.cbegin(), points.cend(), algo);
	int cpt = 0; //Cpt corresponds to the number of straight lines in path
	double length = 0;
	for (auto er = s.begin(), erE = s.end(); er!=erE; ++er) {
		SegmentComputer current(*er);
		Z3i::Point direction;
		PointVector<3, double> intercept, thickness;
		current.getParameters(direction, intercept, thickness);
		if (direction != Z3i::Point::zero) {
			cpt++;
			length += euclideanDistance(*(current.end()-1), *(current.begin()));
		}
	}
	return length;
}

/////////////////////////////////////////////////////////////////////////
//////////////// MAIN ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (skeleton)")
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
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	
	
	
	typedef Z3i::Space Space;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef functors::IntervalForegroundPredicate<Image> Binarizer;
	
	typedef Z3i::KSpace KSpace;
	typedef DGtal::PointVector<3, double> Vector3D;
	typedef DGtal::PointVector<3, int> Point3D;
	typedef LightImplicitDigitalSurface<KSpace, Z3i::DigitalSet >   MyDigitalSurfaceContainer;
	typedef SurfacePoint<Z3i::Point, SCell> SurfacePoint;
	typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
	typedef WeightedDigitalSurface<MyDigitalSurfaceContainer, SurfacePoint> MyWeightedDigitalSurface;
	typedef RadialSurface<MyDigitalSurfaceContainer> MyRadialSurface;
	typedef MetricAdjacency<Space, 3> Graph;
	typedef BreadthFirstVisitor<MyRadialSurface> MyBreadthFirstVisitor;
	typedef BreadthFirstVisitor<MyWeightedDigitalSurface> MyWeightedBreadthFirstVisitor;
	typedef BreadthFirstVisitor<Graph, std::set<Z3i::Point> > PointBreadthFirstVisitor;
	typedef MyBreadthFirstVisitor::Node MyNode;
	typedef MyBreadthFirstVisitor::Size MySize;

	typedef RosenProffittLengthEstimator<set<SCell>> LengthEstimator;
	typedef MSTTangent<Point3D> Tangent;
	typedef Pencil<Point3D, Tangent, Vector3D> Pencil;
	
	trace.beginBlock("Reading file...");
	Image image = VolReader<Image>::importVol(inputFilename);
	trace.endBlock();
	KSpace ks;
	ks.init( image.domain().lowerBound(), 
			 image.domain().upperBound(), true );
	KSpace::SCellSet boundary;
	Z3i::DigitalSet set3d (image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image> (set3d, image, 
												  thresholdMin, thresholdMax);
	trace.info() << "Finding a bel" << endl;
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 100000 );
	//bel = Z3i::SCell({465,516,289},false);
	typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
	MySurfelAdjacency surfAdj( true );
	MyDigitalSurfaceContainer* ptrSurfContainer = 
		new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
	MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
	//! [volBreadthFirstTraversal-SetUpDigitalSurface]

    set<SCell> surfelSet;
	map<SCell, Point> surfaceVoxelSet;
	DigitalSet setPredicate(image.domain());
	set<Point> surfacePointSet;
	for (auto it = digSurf.begin(), itE = digSurf.end(); it!=itE; ++it) {
		surfelSet.insert(*it);
		setPredicate.insert(ks.sCoords(*it));
	}
	vector<SurfacePoint> weightedSurfaceV = SurfaceUtils::computeSurfelWeight<SurfacePoint>(ks, set3d, surfelSet, surfacePointSet, surfaceVoxelSet);
	MyWeightedDigitalSurface weightedSurface(digSurf, weightedSurfaceV);
	MyNode node;
	bool isPathFound = false;
	
    trace.beginBlock("Compute path");
	Z3i::SCell end;
	int numberToFind = 1;
	int i = 0;
	vector<Pencil> pencils;
	set<SCell> thePath;
		
	QApplication application(argc,argv);
	Viewer3D<> viewer( ks);
	viewer.show();
	while (i < numberToFind) {
		bel = Surfaces<KSpace>::findABel( ks, set3d, 100000 );
		auto itDirs  = ks.sDirs(bel);
		
		const int notDirection = *itDirs;
		MyRadialSurface radSurf(digSurf, notDirection);
		auto tracker = digSurf.container().newTracker( bel );
		tracker->move(bel);
		SCell s;
		
//		bel = {{59,32,131}, false};
//		bel = {{27,10,115}, false};
//	  	bel = {{109,80,189}, false};
		viewer << CustomColors3D(Color::Green, Color::Green) << bel;
		MyBreadthFirstVisitor visitor( radSurf, bel );
		isPathFound = false;
		i++;
		map<Z3i::SCell, Z3i::SCell> aMapPrevious;
		trace.info() << bel << endl;
		int stepForShortestPath = -2;
		double distanceForShortestPath = numeric_limits<double>::max();
		Vector3D theNormal;
		
		map<SCell, int> mapVisit;
		
		while (!visitor.finished() && !isPathFound) {
			mapVisit[node.first] = node.second;
			node = visitor.current();
			auto itDirs  = ks.sDirs(node.first);
			tracker->move(node.first);
			vector<Z3i::SCell> neighbors = computeMapInOneDirection(s, tracker, itDirs, notDirection);

			for (auto it = neighbors.begin(), itE = neighbors.end(); it != itE; ++it) {
				if (*it == bel) trace.info() << node.second << endl;
				// auto itMapExisting = aMapPrevious.find(*it);
				// if (itMapExisting == aMapPrevious.end())
				// {
				// 	aMapPrevious[*it] = node.first;
				// }
				// else {
				// 	Z3i::SCell tmp = node.first;
				// 	Z3i::SCell tmp2 = aMapPrevious[*it];
				// 	vector<Z3i::Point> aDirectionV;
				// 	vector<Z3i::Point> anotherDirectionV;
				// 	set<Z3i::SCell> path1;
				//     set<Z3i::SCell> path2;
				
				// 	while (tmp != bel) {
				// 		path1.insert(tmp);
				// 		Z3i::SCell previous = aMapPrevious[tmp];
				// 		aDirectionV.push_back(ks.sCoords(tmp) - ks.sCoords(previous));
				// 		tmp = previous;
				// 	}
			
				// 	while (tmp2 != bel) {
				// 		path2.insert(tmp2);
				// 		Z3i::SCell previous = aMapPrevious[tmp2];
				// 		anotherDirectionV.push_back(ks.sCoords(tmp2) - ks.sCoords(previous));
				// 		tmp2 = previous;
				// 	}
				//     set<Z3i::SCell> x;
				// 	set_union(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(x, x.end()));
							  
				// 	Z3i::Point center = (ks.sCoords(*it) + ks.sCoords(bel)) / 2;
				// 	//check both paths going back			
				// 	vector<Z3i::SCell> intersection;
				// 	for (auto it = path1.begin(), itE = path1.end(); it!=itE; ++it) {
				// 		for (auto it2 = path2.begin(), itE2 = path2.end(); it2!=itE2; ++it2) {
				// 			if (*it == *it2) {
				// 				intersection.push_back(*it);
				// 			}
				// 		}
				// 	}
				// 	if (intersection.size() == 0) {
				// 	}
				// }
			}
//			viewer << CustomColors3D(Color::Green, Color::Green) << node.first;
			visitor.expand();
		}
		for (auto it = thePath.begin(), itE = thePath.end(); it != itE; ++it) {
			viewer << CustomColors3D(Color::Red, Color::Red) << *it;
		}

		typedef typename DGtal::HueShadeColorMap<double, 2 > ColorMap;
		const ColorMap colormap(1, node.second);
		
		auto visitedV = visitor.visitedVertices();
		for (auto it = mapVisit.begin(), itE = mapVisit.end();
			 it != itE; ++it) {
			viewer << CustomColors3D(colormap(it->second), colormap(it->second)) << it->first;
		}
	}
	trace.endBlock();
	
	
	//SliceUtils::slicesFromPlanes(viewer, pencils, image, "img/slice");
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	for (auto it = set3d.begin(); it != set3d.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
	}  
	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();

	return 0;
}
