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
 * @author Jacques-Olivier Lachaud (\c jacques-ol// ivier.lachaud@univ-savoie.fr )
//  * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
//  *
//  * @date 2014/01/31
//  *
//  * Computes the 2d voronoi map of a list of digital points.
//  *
//  * This file is part of the DGtal library.
//  */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/kernel/sets/DigitalSetBySTLSet.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/geometry/curves/estimation/RosenProffittLocalLengthEstimator.h"
#include "DGtal/geometry/curves/estimation/DSSLengthEstimator.h"
#include "DGtal/geometry/curves/estimation/BLUELocalLengthEstimator.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/geometry/curves/estimation/FPLengthEstimator.h"


///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

namespace po = boost::program_options;


typedef Z3i::Space Space;
typedef Z3i::KSpace KSpace;
typedef Z3i::Point Point;
typedef Z3i::RealPoint RealPoint;
typedef Z3i::RealVector RealVector;
typedef HyperRectDomain<Space> Domain;
typedef KSpace::Surfel Surfel;
typedef KSpace::Cell Cell;

typedef ImageContainerBySTLMap<Domain, int> Image;
typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;



typedef Z3i::Curve::ArrowsRange Range;
typedef vector<pair<Point, Z3i::Vector> >::const_iterator ArrowIterator;
typedef vector<Point>::const_iterator PointIterator;
typedef RosenProffittLocalLengthEstimator<ArrowIterator> RosenLength;
typedef DSSLengthEstimator<PointIterator> DSSLength;
typedef BLUELocalLengthEstimator<ArrowIterator> BLUELengthEstimator;
typedef  DigitalSetBySTLSet<DGtal::HyperRectDomain<DGtal::SpaceND<3u, int> >, std::less<DGtal::PointVector<3u, int, boost::array<int, 3ul> > > >::Iterator Iterator;
typedef FPLengthEstimator<PointIterator> FPLocalLengthEstimator;
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
	const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
//QApplication application(argc,argv);
	

// parse command line ----------------------------------------------
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
		("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" ) 
		("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" );    

	bool parseOK=true;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, general_opt), vm);  
	}catch(const std::exception& ex){
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
  
	if(! vm.count("input"))
	{
		trace.error() << " The file name was defined" << endl;      
		return 0;
	}
	string inputFilename2 = vm["input"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
  
	string inputFilename = examplesPath + inputFilename2;
	trace.info() << "File             = " << inputFilename << std::endl;
	trace.info() << "Min image thres. = " << thresholdMin << std::endl;
	trace.info() << "Max image thres. = " << thresholdMax << std::endl;

	KSpace ks;

// Reads the volume
	trace.beginBlock( "Loading image into memory and build digital surface." );
	Image image = GenericReader<Image>::import(inputFilename);
	trace.info() << image.domain().upperBound() << endl;
	ks.init(image.domain().lowerBound(), image.domain().upperBound(), false);
	Z3i::DigitalSet set3d (image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image>(set3d, image, thresholdMin, thresholdMax);
	SurfelAdjacency<KSpace::dimension> SAdj( true );
	vector<vector<Z3i::SCell> > vectorSCell;
	vector<pair<Point, Z3i::Vector> > vectorArrows;
	vector<Point> points;

	/** Filling point vector **/
	for (Iterator iterator = set3d.begin(); iterator != set3d.end(); ++iterator) {
		points.push_back((*iterator));
	}
	
     // Getting the consecutive surfels of the 2D boundary
	Surfaces<KSpace>::extractAllConnectedSCell(vectorSCell, ks, SAdj, set3d);
	trace.endBlock();
	trace.info() << vectorSCell.size() << endl;

	KSpace ks2;
	ks2.init(Point(-25,-25,-25), Point(25,25,25), false);
	functors::SCellToArrow<KSpace> convertS2A(ks2);
	for (vector<vector<Z3i::SCell> >::iterator i = vectorSCell.begin();
		 i != vectorSCell.end(); ++i) {
		for (vector<Z3i::SCell>::iterator it = (*i).begin();
			 it != (*i).end(); ++it) {
			vectorArrows.push_back(convertS2A((*it)));
		}
	}
	
	//adjacency
//! [imageGridCurveEstimator-prepareTracking]

//! [imageGridCurveEstimator-tracking]
//extraction of all the contours;
    ArrowIterator it = vectorArrows.begin();
	ArrowIterator itE = vectorArrows.end();
//trace.info() << (*it).first << endl;

	/*DSSLength dssLength;
	dssLength.init(1, points.begin(), points.end());
	trace.info() << "DSS length= " << dssLength.eval() << endl;*/


	FPLocalLengthEstimator fpLength;
	fpLength.init(1, points.begin(), points.end(), false);
	trace.info() << "FP longueur= " << fpLength.eval() << endl;

	BLUELengthEstimator blueLength;
	blueLength.init(1, vectorArrows.begin(), vectorArrows.end(), false);
	trace.info() << "BLUE longueur= " << blueLength.eval() << endl;

    RosenLength rosenLength;
    rosenLength.init(1, it, itE, false);
//	length.eval();
	trace.info() << "Rosen Longueur= " << rosenLength.eval() << endl;
    
}
