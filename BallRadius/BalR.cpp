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
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/topology/CanonicSCellEmbedder.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"

using namespace DGtal;
using namespace std;
using namespace Z3i;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GreedySegmentation
///////////////////////////////////////////////////////////////////////////////

template <typename SCellSet, typename SurfelPredicate >
void constructSubGraph( SCellSet & surface,
						const KSpace & K,
						const SurfelAdjacency<KSpace::dimension> & surfel_adj,
						const SurfelPredicate & sp,
						const SCell & start_surfel,
						const SCell & end_surfel,
						Viewer3D<>& viewer)
{
	BOOST_CONCEPT_ASSERT(( concepts::CSurfelPredicate<SurfelPredicate> ));

	SCell b;  // current surfel
	SCell bn; // neighboring surfel
	ASSERT( K.sIsSurfel( start_surfel ) );
	surface.clear(); // boundary being extracted.

	SurfelNeighborhood<KSpace> SN;
	SN.init( &K, &surfel_adj, start_surfel );
	std::queue<SCell> qbels;
	qbels.push( start_surfel );
	surface.insert( start_surfel );
	// For all pending bels
	map<Z3i::SCell, Z3i::SCell> aMapPrevious;
	while ( ! qbels.empty() )
    {
		b = qbels.front();
		qbels.pop();
		SN.setSurfel( b );
		if (b == end_surfel) {
			vector<SCell> path;
			SCell previous = end_surfel;
			while (previous != start_surfel) {
				path.push_back(previous);
				viewer << CustomColors3D(Color::Black, Color::Black) << previous;
				previous = aMapPrevious[previous];
			}
			cout << "Other path size" << path.size() << endl;
		}
		for ( auto q = K.sDirs( b ); q != 0; ++q )
        {
			Dimension track_dir = *q;
			// ----- 1st pass with positive orientation ------
			if ( SN.getAdjacentOnSurfelPredicate( bn, sp, track_dir, true ) )
            {
				if ( surface.find( bn ) == surface.end() )
                {
					surface.insert( bn );
					qbels.push( bn );
					aMapPrevious[bn] = b;
                }
            }
			// ----- 2nd pass with negative orientation ------
			if ( SN.getAdjacentOnSurfelPredicate( bn, sp, track_dir, false ) )
            {
				if ( surface.find( bn ) == surface.end() )
                {
					surface.insert( bn );
					qbels.push( bn );
					aMapPrevious[bn] = b;
                }
            }
        } // for ( DirIterator q = K.sDirs( b ); q != 0; ++q )
    } // while ( ! qbels.empty() )
}


float deduceCoordinateFromEquationThreeUnknown(float a, float b, float c, float first_value, float second_value) {
	return (-(a * first_value + b * second_value) / c); 
}




template <typename Pencil, typename Iterator>
vector<Pencil> orthogonalPlanesWithNaiveTangents(Iterator itb, Iterator ite) {
	typedef typename Pencil::P Vector;
	vector<Pencil> tangents;
	for (; itb != ite; ++itb) {
		Iterator nextIt = itb;
		++nextIt;
		Vector naiveTangent = (*nextIt) - (*itb);
		tangents.push_back(Pencil(*itb, naiveTangent));
	}
	return tangents;
}

template <typename Point>
bool multipleDirection(const vector<Point>& aVector, vector<Point>& inflexionPoint) {
	bool xOrientation = false, yOrientation = false, zOrientation = false;
	bool posx = false, posy = false, posz = false, negx = false, negy = false, negz = false;
	for (int i = 0; i < aVector.size() - 2; i++) {
		Point diff = aVector[i] - aVector[i+1];
		Point toPush = aVector[i+1];
		trace.info() << aVector[i] << " " << aVector[i+1] << endl;
		if (diff[0] > 0 && !posx) {
			posx = true;
			inflexionPoint.push_back(toPush);
		}
		if (diff[0] < 0 && !negx) {
			negx = true;
			inflexionPoint.push_back(toPush);
		}
		if (diff[1] > 0 && !posy) {
		    posy = true;
			inflexionPoint.push_back(toPush);
		}
		if (diff[1] < 0 && !negy) {
			negy = true;
			inflexionPoint.push_back(toPush);
		}
		if (diff[2] > 0 && !posz) {
			posz = true;
			inflexionPoint.push_back(toPush);
		}
		if (diff[2] < 0 && !negz) {
			negz = true;
			inflexionPoint.push_back(toPush);
		}
	}
	xOrientation = posx && negx;
	yOrientation = posy && negy;
	zOrientation = posz && negz;
	return (xOrientation && yOrientation) || (xOrientation && zOrientation) || (zOrientation && yOrientation);
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
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 1000000 );

	trace.info() << bel << endl;
//	bel = Z3i::SCell({-12,1,107},true);
	typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
	MySurfelAdjacency surfAdj( true );
	MyDigitalSurfaceContainer* ptrSurfContainer = 
		new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
	MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
	//! [volBreadthFirstTraversal-SetUpDigitalSurface]

	KSpace::SCellSet surfelSet;
	set<Point> surfaceVoxelSet;
	DigitalSet setPredicate(image.domain());
	for (auto it = digSurf.begin(), itE = digSurf.end(); it!=itE; ++it) {
		surfelSet.insert(*it);
		surfaceVoxelSet.insert(ks.sCoords(*it));
		setPredicate.insert(ks.sCoords(*it));
	}
	functors::SurfelSetPredicate<KSpace::SCellSet, SCell> setP(surfelSet);
	
	//! [volBreadthFirstTraversal-ExtractingSurface]
	trace.beginBlock( "Extracting boundary by tracking from an initial bel." );
	Graph graph;
	MyBreadthFirstVisitor visitor( digSurf, bel );
	MyNode node;
	 
	bool isPathFound = false;
	map<Z3i::SCell, pair<Z3i::SCell, unsigned int> > aMapPrevious;
	viewer << CustomColors3D(Color::Green, Color::Green) << bel;
	
    trace.beginBlock("Compute path");
	Z3i::SCell end;
	while (!visitor.finished() && !isPathFound) {
		node = visitor.current();
		vector<Z3i::SCell> neighbors;
		back_insert_iterator<vector<Z3i::SCell> > iter(neighbors);
	    visitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), itE = neighbors.end(); it != itE; ++it) {
			
			auto itMapExisting = aMapPrevious.find(*it);
			if (itMapExisting == aMapPrevious.end())
		    {
				aMapPrevious[*it] = pair<Z3i::SCell, unsigned int>(node.first, node.second+1);			
				
			}
			else {
				Z3i::SCell tmp = node.first;
				Z3i::SCell tmp2 = aMapPrevious[*it].first;
				vector<Z3i::Point> aDirectionV;
				vector<Z3i::Point> anotherDirectionV;
				set<Z3i::SCell> path1;
				set<Z3i::SCell> path2;
				while (tmp != bel) {
				    path1.insert(tmp);
					Z3i::SCell previous = aMapPrevious[tmp].first;
					aDirectionV.push_back(ks.sCoords(tmp) -ks.sCoords(previous));
					tmp = previous;
					
				}
			
				while (tmp2 != bel) {
					path2.insert(tmp2);
					Z3i::SCell previous = aMapPrevious[tmp2].first;
					anotherDirectionV.push_back(ks.sCoords(tmp2));
					tmp2 = previous;
				}
				float size1 = 0;
				float size2 = 0;
				float increment = 2;
				for (int i = 1; i < aDirectionV.size() - 1; i++) {
					if (aDirectionV[i] == aDirectionV[i-1] && aDirectionV[i] == aDirectionV[i+1]) 
						size1 += increment;
					else
						size1++;
					if (anotherDirectionV[i] == anotherDirectionV[i-1] && anotherDirectionV[i] == anotherDirectionV[i+1])
						size2 += increment;
					else
						size2++;
				}
				
				if (size1 == size2) {
					

					//check shorter path
					vector<SCell> shorterPath;
				
					//check both paths going back
				
				
					set<Z3i::SCell> intersection;
					set_intersection(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(intersection, intersection.begin()));
					vector<Point> inflexDirectionV;
					vector<Point> anotherInflexDirectionV;
					bool isPath1 = multipleDirection(aDirectionV, inflexDirectionV);
					bool isPath2 = multipleDirection(anotherDirectionV, anotherInflexDirectionV);
					// vector<Point> concatenation;
					// concatenation.reserve(aDirectionV.size() + anotherDirectionV.size());
					// concatenation.insert(concatenation.end(), aDirectionV.begin(), aDirectionV.end());
					// concatenation.insert(concatenation.end(), anotherDirectionV.begin(), anotherDirectionV.end());
					// bool isPath = multipleDirection(concatenation, inflexDirectionV);
					if (intersection.size() == 0 && isPath1 && isPath2) {
						
					
						for (auto it = aDirectionV.begin(), itE = aDirectionV.end(); it != itE; ++it) {
							trace.info() << *it << endl;
						}
							   
						viewer << CustomColors3D(Color::Yellow, Color::Yellow) << *path1.lower_bound(bel) << *path2.lower_bound(bel);
						end = *it;
						KSpace::SCellSet aSCellSet;
						constructSubGraph(aSCellSet, ks, surfAdj, setP, bel, end, viewer);
						cout << shorterPath.size() << endl;
						for (auto it = shorterPath.begin(), itE = shorterPath.end(); it != itE; ++it) {
							viewer << CustomColors3D(Color::Black, Color::Black) << *it;
						}
						trace.info() << "path1 size " << path1.size() << " path2 size " << path2.size() << "distance " << node.second << endl;
						isPathFound = true;
						Z3i::SCell start = *it;
						viewer << CustomColors3D(Color::Green, Color::Green) << node.first;
						viewer << CustomColors3D(Color::Yellow, Color::Yellow) << start;
						while (start != bel) {
							start  = aMapPrevious[start].first;
							viewer << CustomColors3D(Color::Red, Color::Red) << start;
						}
						Z3i::SCell otherStart = node.first;
						trace.info() << otherStart << *it << endl;
						while (otherStart != bel) {
							viewer << CustomColors3D(Color::Aqua, Color::Aqua) << otherStart;
							otherStart = aMapPrevious[otherStart].first;
					
						}
						break;
					}
				}
				
			}
			
		}
		
		visitor.expand();
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
