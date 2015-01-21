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
// DSS LENGTH
///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <QtGui/qapplication.h>
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
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/geometry/curves/FP.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "StandardDSS26Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

namespace po = boost::program_options;


typedef Z3i::Space Space;
typedef Z3i::KSpace KSpace;
typedef KhalimskySpaceND<3, int> K3;
typedef Z3i::Point Point;
typedef Z3i::RealPoint RealPoint;
typedef Z3i::RealVector RealVector;
typedef HyperRectDomain<Space> Domain;
typedef KSpace::Surfel Surfel;
typedef KSpace::Cell Cell;
typedef ImageContainerBySTLMap<Domain, int> Image;
typedef vector<Point>::const_iterator PointIterator;
typedef  DigitalSetBySTLSet<DGtal::HyperRectDomain<DGtal::SpaceND<3u, int> >, std::less<DGtal::PointVector<3u, int, boost::array<int, 3ul> > > >::Iterator Iterator;


const Color  AXIS_COLOR( 0, 0, 0, 255 );
const double AXIS_LINESIZE = 0.1;
const Color  XY_COLOR( 0, 0, 255, 50 );
const Color  XZ_COLOR( 0, 255, 0, 50 );
const Color  YZ_COLOR( 255, 0, 0, 50 );


template <typename KSPace, typename Space>
double estimateLength(const vector<Point> & points, Viewer3D<Space, KSPace> & viewer, bool display=false) {
	typedef typename PointIterator::value_type Point3d;
	typedef StandardDSS6Computer<PointIterator,int,8> SegmentComputer;
	typedef SaturatedSegmentation<SegmentComputer> Decomposition;
	typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
	HueShadeColorMap<int> cmap_hue( 0, 3, 1 );
	SegmentComputer algo;
	Decomposition decomp(points.begin(), points.end(), algo);
	double totalEstimatedLength= 0.0;
	viewer << SetMode3D( algo.className(), "BoundingBox" );
	int c = 0;
	for (SegmentComputerIterator i = decomp.begin(); i != decomp.end(); ++i) {
		SegmentComputerIterator tmp = i;
		SegmentComputer ms3d(*i);
		Point3d first = *ms3d.begin();
		Point3d last = *(ms3d.end() -1);
		trace.info() << "-  MS3D,"
					 << " [" << first[ 0 ] << "," << first[ 1 ] << ","<< first[ 2 ] << "]"
					 << "->[" << last[ 0 ] << "," << last[ 1 ] << ","<< last[ 2 ] << "]" << endl;
		double distanceDSS = sqrt(pow((last[0] - first[0]), 2) + pow((last[1] - first[1]), 2) + pow((last[2] - first[2]), 2));
		totalEstimatedLength += distanceDSS;
		Color color3d = cmap_hue(c);
		if (display) {
			viewer << CustomColors3D(Color::Red, Color::Red) << ms3d;
		}
		c++;
	}
	return totalEstimatedLength;
}



template <typename KSpace, typename StandardDSS6Computer, typename space, typename kspace >
void displayProj2d( Viewer3D<space, kspace> & viewer,
		    const KSpace & ks, const StandardDSS6Computer & dss3d,
		    const DGtal::Color & color2d )
{
  typedef typename StandardDSS6Computer::ArithmeticalDSSComputer2d ArithmeticalDSSComputer2d;
  typedef typename ArithmeticalDSSComputer2d::ConstIterator ConstIterator2d;
  typedef typename ArithmeticalDSSComputer2d::Point Point2d;
  typedef typename KSpace::Cell Cell;
  typedef typename KSpace::Point Point3d;
  Point3d b = ks.lowerBound();
  for ( DGtal::Dimension i = 0; i < 3; ++i )
    {
      const ArithmeticalDSSComputer2d & dss2d = dss3d.arithmeticalDSS2d( i );
      for ( ConstIterator2d itP = dss2d.begin(), itPEnd = dss2d.end(); itP != itPEnd; ++itP )
	{
	  Point2d p = *itP;
	  Point3d q;
	  switch (i) {
	  case 0: q = Point3d( 2*b[ i ]  , 2*p[ 0 ]+1, 2*p[ 1 ]+1 ); break;
	  case 1: q = Point3d( 2*p[ 0 ]+1, 2*b[ i ]  , 2*p[ 1 ]+1 ); break;
	  case 2: q = Point3d( 2*p[ 0 ]+1, 2*p[ 1 ]+1, 2*b[ i ]   ); break;
	  }
	  Cell c = ks.uCell( q );
	  trace.info() << c << endl;
//	  viewer << CustomColors3D( color2d, color2d ) << c;
	}
    }
}

/*template <typename KSPace, typename Space>
double estimateLengthWithFP(const vector<Point> &points, Viewer3D<Space, KSPace> & viewer, bool display = false) {
	typedef FP<PointIterator, int, 4> MariannePolygon;
	double totalEstimatedLength = 0.0;
	MariannePolygon polygonalisation(points.begin(), points.end());
	for (list<Z2i::Point>::const_iterator it = polygonalisation.polygon().begin(); it != polygonalisation.polygon().end(); ++it) {
		list<Z2i::Point>::const_iterator next = it;
		++next;
		if (next != polygonalisation.polygon().end()) {
			Z2i::Point first = *it;
			Z2i::Point last = *next;
			trace.info() << "-  FP,"
						 << " [" << first[ 0 ] << "," << first[ 1 ] << ","<< first[ 2 ] << "]"
						 << "->[" << last[ 0 ] << "," << last[ 1 ] << ","<< last[ 2 ] << "]" << endl;
			double distanceFP = sqrt(pow((last[0] - first[0]), 2) + pow((last[1] - first[1]), 2) + pow((last[2] - first[2]), 2)); 
			totalEstimatedLength += distanceFP;
		}
	}
	return totalEstimatedLength;
	}*/



///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
	const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
	QApplication application(argc,argv);
	Viewer3D<> viewer;
	

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


// Reads the volume
	trace.beginBlock( "Loading image into memory and build digital surface." );
	vector<Point> points=  PointListReader<Point>::getPointsFromFile(inputFilename); 
	/*Image image = VolReader<Image>::importVol(inputFilename);
//trace.info() << image.domain().upperBound() << endl;
ks.init(image.domain().lowerBound(), image.domain().upperBound(), false);
//Z3i::DigitalSet set3d (image.domain());
//SetFromImage<Z3i::DigitalSet>::append<Image>(set3d, image, thresholdMin, thresholdMax);
*/


    K3 ks;
	Point lowerBound = points[ 0 ];
	Point upperBound = points[ 0 ];
	for ( unsigned int j = 1; j < points.size(); ++j )
    {
		lowerBound = lowerBound.inf( points[ j ] );
		upperBound = upperBound.sup( points[ j ] );
    }
	lowerBound -= Point::diagonal( 3 );
	upperBound += Point::diagonal( 4 );
	ks.init(lowerBound, upperBound, true);
		
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	viewer.show();
	for (vector<Point>::iterator it = points.begin(); it != points.end(); ++it) {
		trace.info() << (*it) << endl;
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
	}
	const string mode = "WIRED";
	//displayAxes<Point, SpaceND<3, int>::RealPoint, Z3i::Space, Z3i::KSpace>(viewer, ks.lowerBound(), ks.upperBound(), mode);
	trace.info() << "DSS estimated length = " << estimateLength<>(points, viewer, true) << endl;

	viewer << Viewer3D<>::updateDisplay;
	application.exec();
}
