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
#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "geometry/VoronoiCovarianceMeasure.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/graph/DepthFirstVisitor.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/TangentUtils.h"
#include "geometry/SliceUtils.h"
#include "geometry/Pencil.h"
#include "geometry/MSTTangent.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
int main( int  argc, char**  argv )
{
	QApplication application(argc,argv);

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
	typedef functors::HatPointFunction<Point,double> KernelFunction;
  
	typedef MSTTangent<Point> Tangent;
	typedef Pencil<Point, Tangent, RealPoint> Pencil;
	typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("skeleton,s", po::value<std::string>(), "vol file (skeleton)")
		("input,i", po::value<std::string>(), "vol file (corresponding volume)")
		("output,o", po::value<std::string>(), "sliced vol file with orthogonal planes")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		("radiusInside,r", po::value<double>()->default_value(10), "radius of the ball inside voronoi cell")
		("radiusNeighbour,R", po::value<double>()->default_value(10), "radius of the ball for the neighbourhood")
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
	string skeletonFilename = vm["skeleton"].as<std::string>();
	string inputFilename = vm["input"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	double R = vm["radiusInside"].as<double>();
	double r = vm["radiusNeighbour"].as<double>();


	Image volume = VolReader<Image>::importVol(inputFilename);
	Z3i::Domain domainVolume = volume.domain();
  
	Image image = VolReader<Image>::importVol(skeletonFilename); 

	Point p;	
	vector<Point> vPoints;
	for (auto it = image.domain().begin(), itE = image.domain().end(); it != itE; ++it) {
		if (image(*it) >= thresholdMin && image(*it) <= thresholdMax) {
			p = *it;
		}
	}
// We have to visit the direct neighbours in order to have a container with voxels
// ordered sequentially by their connexity
// Otherwise we have a point container with points which are not neighbours
// and this impairs maximal segment recognition
	typedef MetricAdjacency<Space, 3> Graph;
	typedef DepthFirstVisitor<Graph, set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;
  
	Graph graph;
	Visitor visitor( graph, p );
	MyNode node;

    
	while ( !visitor.finished() ) 
	{
		node = visitor.current();
		if ( image.domain().isInside(node.first) &&
			 image(node.first) >= thresholdMin &&
			 image(node.first) <= thresholdMax ) { //is inside domain
			vPoints.push_back(node.first);
			visitor.expand();
		}
		else
			visitor.ignore();
	}

	
	trace.info() << "Big radius   R = " << R << std::endl;
	trace.info() << "Small radius r = " << r << std::endl;
  
	const double size = 20.0; // size of displayed normals


//Computing lambda MST tangents
	std::vector<Pencil> tangents = TangentUtils::orthogonalPlanesWithTangents<Pencil>(vPoints.begin(), vPoints.end());
  
	Metric l2;
	VCM vcm( R, ceil( r ), l2, true );
	vcm.init( vPoints.begin(), vPoints.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, r );

// Flat zones are metallic blue, slightly curved zones are white,
// more curved zones are yellow till red.
	Matrix vcm_r, evec;
	RealVector eval;
	Viewer3D<> viewer;
	viewer.show();
	int sliceNumber = 0;
 
	const int IMAGE_PATCH_WIDTH = 100;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH), volume.domain().upperBound() + Z3i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
									  DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::functors::Identity idV;
	
	for ( auto it = ++tangents.begin(), itE = tangents.end();
		  it != itE; ++it )
	{
		DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, it->getPoint(), it->getTangent(), IMAGE_PATCH_WIDTH);
		ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
		extractedImage.setDefaultValue(0);
		double radius  = SliceUtils::computeRadiusFromImage(extractedImage, thresholdMin, thresholdMax);
		// Compute VCM and diagonalize it.
		viewer.setFillColor(Color::Gray);
		viewer.setFillTransparency(255);
		
		viewer << it->getPoint();
		/*if (radius > 0) {
			vcm.updateProximityStructure(radius, vPoints.begin(), vPoints.end());
			chi = KernelFunction( 1.0, radius);
			}*/

		vcm_r = vcm.measure( chi, it->getPoint() );
		LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
		double feature = eval[1] /  (eval[0] + eval[ 1 ]+  eval[3]);
	    // Display normal
		RealVector n = evec.column( 2 );
		n*=size;
		RealPoint p( it->getPoint()[ 0 ], it->getPoint()[ 1 ], it->getPoint()[2] );
		RealVector n2 = evec.column( 1 );
		n2*=size;

		if(sliceNumber > 2280 && sliceNumber < 2420 && sliceNumber % 10 ==0) {
			viewer.setLineColor(Color::Blue);
			viewer.setFillColor(Color::Blue);
			viewer.setFillTransparency(150);
			viewer.addQuad(p-n-n2,p-n+n2,p+n+n2,p+n-n2);
		}
		sliceNumber++;
	}

	for (auto it = volume.domain().begin(), ite = volume.domain().end(); it != ite; ++it) {
		if (volume(*it) >= thresholdMin)
			viewer << CustomColors3D(Color(255,0,0,20), Color(255,0,0,20))<<*it;
	}
  
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////