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
 * @date 2014/01/31
 *
 * Computes the 2d voronoi map of a list of digital points.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <QtGui/qapplication.h>

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
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/GenericReader.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;
///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
	QApplication application(argc,argv);

	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef Z3i::Point Point;
	typedef Z3i::RealPoint RealPoint;
	typedef Z3i::RealVector RealVector;
	typedef HyperRectDomain<Space> Domain;
	typedef KSpace::Surfel Surfel;
	typedef KSpace::Cell Cell;

	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef ImplicitDigitalSurface< KSpace, ThresholdedImage > DigitalSurfaceContainer;

	//! [DVCM3D-typedefs]
	typedef ExactPredicateLpSeparableMetric<Space, 2> Metric;          // L2-metric type
	typedef functors::HatPointFunction<Point,double>  KernelFunction;  // chi function type 
	typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric,
													  KernelFunction > VCMOnSurface;
	typedef VCMOnSurface::Surfel2Normals::const_iterator S2NConstIterator;
	//! [DVCM3D-typedefs]
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
		("featureThreshold,T", po::value<double>()->default_value(0.1), "Feature for sharp edges threshold")
		("thresholdMin,m",  po::value<int>()->default_value(0), "threshold min to define binary shape" ) 
		("thresholdMax,M",  po::value<int>()->default_value(255), "threshold max to define binary shape" )    
		("transparency,t",  po::value<uint>()->default_value(255), "transparency") ; 

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
	if(!vm.count("input"))
	{
		trace.error() << " The file name was defined" << endl;      
		return 0;
	}
	string inputFilename = vm["input"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	unsigned char transp = vm["transparency"].as<uint>();
	double T = vm["featureThreshold"].as<double>();
 
	trace.info() << "File             = " << inputFilename << std::endl;
	trace.info() << "Min image thres. = " << thresholdMin << std::endl;
	trace.info() << "Max image thres. = " << thresholdMax << std::endl;
	const double R = 5;
	trace.info() << "Big radius     R = " << R << std::endl;
	const double r = 3;
	trace.info() << "Small radius   r = " << r << std::endl;
	const double trivial_r = 3;
	trace.info() << "Trivial radius t = " << trivial_r << std::endl; // for orienting the directions given by the tensor.
	trace.info() << "Feature thres. T = " << T << std::endl; // threshold for displaying features as red.
	

	const double size = 1.0; // size of displayed normals.

	KSpace ks;
	// Reads the volume
	trace.beginBlock( "Loading image into memory and build digital surface." );
	Image image = GenericReader<Image>::import(inputFilename );
	ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
	trace.endBlock();
	trace.beginBlock( "Extracting boundary by scanning the space. " );
	ks.init( image.domain().lowerBound(),
			 image.domain().upperBound(), true );
	SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
	Surfel bel = Surfaces<KSpace>::findABel( ks, thresholdedImage, 10000 );
	DigitalSurfaceContainer* container = 
		new DigitalSurfaceContainer( ks, thresholdedImage, surfAdj, bel, false  );
	DigitalSurface< DigitalSurfaceContainer > surface( container ); //acquired
	trace.info() << "Digital surface has " << surface.size() << " surfels." << std::endl;
	trace.endBlock();

	//! [DVCM3D-instantiation]
	Surfel2PointEmbedding embType = Pointels; // Could be Pointels|InnerSpel|OuterSpel; 
	Metric l2;                                // Euclidean L2 metric 
	KernelFunction chi( 1.0, r );             // hat function with support of radius r
	VCMOnSurface vcm_surface( surface, embType, R, r, 
							  chi, trivial_r, l2, true /* verbose */ );
	//! [DVCM3D-instantiation]

	trace.beginBlock( "Displaying VCM" );
	Viewer3D<> viewer( ks );
	Cell dummy;
	viewer.setWindowTitle("3D VCM viewer");
	viewer << SetMode3D( dummy.className(), "Basic" );
	viewer.show();

	GradientColorMap<double> grad( 0, T );
	grad.addColor( Color( 128, 128, 255 ) );
	grad.addColor( Color( 255, 255, 255 ) );
	grad.addColor( Color( 255, 255, 0 ) );
	grad.addColor( Color( 255, 0, 0 ) );
	RealVector lambda; // eigenvalues of chi-vcm
	for ( S2NConstIterator it = vcm_surface.mapSurfel2Normals().begin(), 
			  itE = vcm_surface.mapSurfel2Normals().end(); it != itE; ++it )
    {
		Surfel s = it->first;
		Point kp = ks.sKCoords( s );
		RealPoint rp( 0.5 * (double) kp[ 0 ], 0.5 * (double) kp[ 1 ], 0.5 * (double) kp[ 2 ] );
		RealVector n = it->second.vcmNormal;
		vcm_surface.getChiVCMEigenvalues( lambda, s );
		double ratio = lambda[ 1 ] / ( lambda[ 0 ] + lambda[ 1 ] + lambda[ 2 ] ); 
		viewer.setFillColor( grad( ratio > T ? T : ratio ) );
		viewer << ks.unsigns( s );
		n *= size;
		viewer.setLineColor( Color::Black );
		viewer.addLine( rp + n, rp - n, 0.1 );
	}
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	trace.endBlock();
	return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
