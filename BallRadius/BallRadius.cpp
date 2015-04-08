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

#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/topology/CanonicSCellEmbedder.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/graph/CVertexPredicate.h"


// Local includes
#include "SurfacePoint.h"
#include "WeightedDigitalSurface.h"
#include "../CreateCurve/Distance.h"
#include "Statistics.h"
using namespace DGtal;
using namespace std;
using namespace Z3i;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GreedySegmentation
///////////////////////////////////////////////////////////////////////////////


template <typename Point>
void visualizePlane(const Point& p, int size, int delta, Viewer3D<> & viewer) {
	int direction = (size - delta) / 2;
	Point n(direction, 0, -delta);
	Point n2(0 , direction, -delta);
	viewer.addQuad(p-n-n2,p-n+n2,p+n+n2,p+n-n2);
}

template <typename Surfel>
void visualizePath(const Surfel& start, const Surfel& end, map<Surfel, Surfel>& previous, Viewer3D<>& viewer, const Color& color) {
	Surfel tmp = end;
	while (tmp != start) {
		viewer << CustomColors3D(color, color) << tmp;
	    tmp = previous[tmp];
	}

}

template <typename Surfel, typename Point, typename KSpace>
bool checkSymmetry(const KSpace& ks, const set<Surfel>& path1, const set<Surfel>& path2, const Point & center) {
	int foundOneSymmetry = 0;
	for (auto it = path1.begin(), itE = path1.end(); it != itE; ++it) {
		Point currentPoint = ks.sCoords(*it);
		Point vectorToCenter = center - currentPoint;
		Point symmetryCurrent = center + vectorToCenter;
		for (auto itS = path2.begin(), itSE = path2.end(); itS != itSE; ++itS) {
			if (symmetryCurrent == ks.sCoords(*itS))
				foundOneSymmetry++;
		}
	}
	return (foundOneSymmetry > (path1.size() / 10) + 1);
}

template <typename Point>
bool multipleDirection(const vector<Point>& aVector, int& cpt) {
	bool xOrientation = false, yOrientation = false, zOrientation = false;
    int posx = 0, posy = 0, posz = 0, negx = 0, negy = 0, negz = 0;
	for (auto it = aVector.begin(), itE = aVector.end(); it != itE; ++it) {
		Point diff = *it;
		if (diff[0] > 0) {
			posx++;
		}
		if (diff[0] < 0) {
			negx++;
		}
		if (diff[1] > 0) {
		    posy++;
		}
		if (diff[1] < 0) {
			negy++;
		}
		if (diff[2] > 0) {
			posz++;
		}
		if (diff[2] < 0) {
			negz++;
		}
	}
	xOrientation = (posx > 0 && negx > 0);
	yOrientation = (posy > 0 && negy > 0);
	zOrientation = (posz > 0 && negz > 0);
	if (xOrientation && yOrientation)
		cpt = posz > 0 ? posz : negz;
	else if (xOrientation && zOrientation)
		cpt = posy > 0 ? posy : negy;
	else if (zOrientation && yOrientation)
		cpt = posx > 0 ? posx : negx;
	return (xOrientation && yOrientation) || (xOrientation && zOrientation) || (zOrientation && yOrientation);
}

template <typename KSpace, typename Surfel, typename Point, typename SurfacePoint>
vector<SurfacePoint> computeSurfelWeight(const KSpace & ks, const DigitalSet & set3d, const set<Surfel>& surfelSet, set<Point> & surfaceVoxelSet) {
	vector<SurfacePoint> weightedSurfaceVector;
	for (auto it = set3d.begin(), itE = set3d.end(); it != itE; ++it) {
		vector<Surfel> aSurfelV;
		SCell current = ks.sSpel(*it);
		int number = 0;
		for (int i = 0; i < 3; i++) {
			auto itSurfel = surfelSet.find(ks.sIncident(current, i, true));
			if (itSurfel != surfelSet.end()) {
				number++;
				aSurfelV.push_back(*itSurfel);
			}
		}

		for (int i = 0; i < 3; i++) {
			auto itSurfel = surfelSet.find(ks.sIncident(current, i, false));
			if (itSurfel != surfelSet.end()) {
				number++;
				aSurfelV.push_back(*itSurfel);
			}
		}
		weightedSurfaceVector.push_back({*it, number, aSurfelV});
		if (number > 0)
			surfaceVoxelSet.insert(ks.sCoords(current));
	}
	return weightedSurfaceVector;
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
	
	typedef PointVector<3,int> Point;
	QApplication application(argc,argv);
	Viewer3D<> viewer;
	viewer.show();
	
	typedef Z3i::Space Space;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef Z3i::KSpace KSpace;
	typedef DGtal::PointVector<3, double> Vector3D;
	typedef DGtal::PointVector<3, int> Point3D;
	typedef LightImplicitDigitalSurface<KSpace, Z3i::DigitalSet >   MyDigitalSurfaceContainer;
	typedef SurfacePoint<Z3i::Point, SCell> SurfacePoint;
	typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
	typedef MetricAdjacency<Space, 2> Graph;
	typedef BreadthFirstVisitor<MyDigitalSurface> MyBreadthFirstVisitor;
	typedef MyBreadthFirstVisitor::Node MyNode;
	typedef MyBreadthFirstVisitor::Size MySize;
	
	
	trace.beginBlock("Reading file...");
	Image image = VolReader<Image>::importVol(inputFilename);
	trace.endBlock();
	KSpace ks;
	ks.init( image.domain().lowerBound(), 
			 image.domain().upperBound(), false );
	KSpace::SCellSet boundary;
	Z3i::DigitalSet set3d (image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image> (set3d, image, 
												  thresholdMin, thresholdMax);
	trace.info() << "Finding a bel" << endl;
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 10000 );
	//bel = Z3i::SCell({465,516,289},false);
	trace.info() << bel << endl;
	typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
	MySurfelAdjacency surfAdj( true );
	MyDigitalSurfaceContainer* ptrSurfContainer = 
		new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
	MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
	//! [volBreadthFirstTraversal-SetUpDigitalSurface]

    set<SCell> surfelSet;
	set<Point> surfaceVoxelSet;
	DigitalSet setPredicate(image.domain());
	
	for (auto it = digSurf.begin(), itE = digSurf.end(); it!=itE; ++it) {
		surfelSet.insert(*it);
		setPredicate.insert(ks.sCoords(*it));
	}
   
	Graph graph;
	
	MyNode node;
	 
	bool isPathFound = false;
	
	viewer << CustomColors3D(Color::Green, Color::Green) << bel;
    trace.beginBlock("Compute path");
	Z3i::SCell end;
	int numberToFind = 10;
	int i = 0;
	while (i < numberToFind) {
		bel = Surfaces<KSpace>::findABel( ks, set3d, 10000 );
		MyBreadthFirstVisitor visitor( digSurf, bel );
		isPathFound = false;
		i++;
		map<Z3i::SCell, Z3i::SCell> aMapPrevious;
		trace.info() << bel << endl;
		while (!visitor.finished() && !isPathFound) {
			if (node.second > 100) isPathFound = true;
			node = visitor.current();
			vector<Z3i::SCell> neighbors;
			back_insert_iterator<vector<Z3i::SCell> > iter(neighbors);
			visitor.graph().writeNeighbors(iter, node.first);
			//viewer << CustomColors3D(Color::White, Color::White) << node.first;
			for (auto it = neighbors.begin(), itE = neighbors.end(); it != itE; ++it) {
			
				auto itMapExisting = aMapPrevious.find(*it);
				if (itMapExisting == aMapPrevious.end())
				{
				
					aMapPrevious[*it] = Z3i::SCell(node.first);			
				
				}
				else {
					Z3i::SCell tmp = node.first;
					Z3i::SCell tmp2 = aMapPrevious[*it];
					vector<Z3i::Point> aDirectionV;
					vector<Z3i::Point> anotherDirectionV;
					set<Z3i::SCell> path1;
					set<Z3i::SCell> path2;
				
					while (tmp != bel) {
						path1.insert(tmp);
						Z3i::SCell previous = aMapPrevious[tmp];
						aDirectionV.push_back(ks.sCoords(tmp) - ks.sCoords(previous));
						tmp = previous;
					}
			
					while (tmp2 != bel) {
						path2.insert(tmp2);
						Z3i::SCell previous = aMapPrevious[tmp2];
						anotherDirectionV.push_back(ks.sCoords(tmp2) - ks.sCoords(previous));
						tmp2 = previous;
					}
			
					/*
					  double size1 = 0;
					  double size2 = 0;
					  float increment = 2.0;
					  for (int i = 1; i < (int)aDirectionV.size() - 1; i++) {
					  if (aDirectionV[i] == aDirectionV[i-1] && aDirectionV[i] == aDirectionV[i+1]) 
					  size1 += increment;
					  else
					  size1++;
					  if (anotherDirectionV[i] == anotherDirectionV[i-1] && anotherDirectionV[i] == anotherDirectionV[i+1])
					  size2 += increment;
					  else
					  size2++;
					  }*/
					//if (size1 == size2) {

					//Checking distance to center
					vector<double> aVector;
					Z3i::Point center = (ks.sCoords(*it) + ks.sCoords(bel)) / 2;
				
					for (auto it = path1.begin(), itE = path1.end(); it != itE; ++it) {
						aVector.push_back(euclideanDistance(ks.sCoords(*it), center));
					}
					double ratio = Statistics::stddev(aVector) / Statistics::mean(aVector);
					//check both paths going back			
					set<Z3i::SCell> intersection;
					
					set_intersection(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(intersection, intersection.begin()));
					int cpt = 0;
					int cpt2 = 0;
					//bool isPath1 = multipleDirection(aDirectionV, cpt);
					//bool isPath2 = multipleDirection(anotherDirectionV, cpt2);
					
					if (ratio < 0.10 && intersection.size() == 0) {
						bool symmetry = checkSymmetry(ks, path1, path2, center);
						if (symmetry){
							viewer << CustomColors3D(Color::Red, Color::Red) << center;
							//All conditions are met: a path is found
							isPathFound = true;

							//Visualize both paths
							viewer << CustomColors3D(Color::Yellow, Color::Yellow) << *it;
							visualizePath(bel, *it, aMapPrevious, viewer, Color::Red);
							visualizePath(bel, node.first, aMapPrevious, viewer, Color::Red);
						
							//Visualize corresponding planes
							/*if (isPath1) {
							  visualizePlane(ks.sCoords(bel), path1.size(), cpt, viewer);
							  } else if (isPath2) {
							  visualizePlane(ks.sCoords(bel), path2.size(), cpt2, viewer);
							  }
							  break;*/
						}
					}
					//}
				
				}
			
			}
			visitor.expand();
		}
	}
	trace.endBlock();
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	for (auto it = set3d.begin(); it != set3d.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
	}  
	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();

	return 0;
}
