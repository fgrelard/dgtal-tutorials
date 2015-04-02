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
#include "DGtal/math/linalg/EigenDecomposition.h"
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
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "../CreateCurve/Ball.h"
#include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/dec/DiscreteExteriorCalculus.h"
#include "DGtal/dec/DiscreteExteriorCalculusSolver.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template <typename Set>
void DEC(const Set& set, const Z3i::Domain & domain, Viewer3D<>& viewer) {
	typedef DiscreteExteriorCalculus<3, EigenLinearAlgebraBackend> Calculus;
	Calculus calculus(set);
	Calculus::PrimalForm0 primal_zero_form(calculus);
	const Calculus::PrimalHodge1 primal_one_hodge = calculus.primalHodge<1>();
	const Calculus::PrimalHodge2 primal_two_hodge = calculus.primalHodge<2>();
	const Calculus::PrimalDerivative1 primal_one_derivative = calculus.derivative<1,PRIMAL>();
	const Calculus::PrimalDerivative0 primal_zero_derivative = calculus.derivative<0, PRIMAL>();
	const Calculus::PrimalForm1 primal_one_form = primal_zero_derivative * primal_zero_form;
	const Calculus::DualForm2 dual_two_form = primal_one_hodge * primal_one_form ;
	//const Calculus::DualVectorField dual_vector_field = calculus.sharp(dual_one_form);
	//const Calculus::PrimalVectorField primal_vector_field = calculus.sharp(primal_one_form);
	//viewer << primal_vector_field;
	Display3DFactory<>::draw(viewer, dual_two_form);
	trace.info() << dual_two_form.getSCell(0) << endl;
/*	Calculus::PrimalVectorField input_vector_field(calculus);
	Display3DFactory<>::draw(viewer, calculus);
	 for (Calculus::Index ii=0; ii<calculus.kFormLength(0, PRIMAL); ii++)
    {
        const Z3i::RealPoint cell_center = Z3i::RealPoint(calculus.getSCell(0, PRIMAL, ii).myCoordinates)/2.;
        input_vector_field.myCoordinates(ii, 0) = 1;
        input_vector_field.myCoordinates(ii, 1) = 0;
        input_vector_field.myCoordinates(ii, 2) = 0;
    }

	 // Display3DFactory<>::draw(viewer, input_vector_field);
	 const Calculus::PrimalForm1 input_one_form = calculus.flat(input_vector_field);
	
	const Calculus::PrimalDerivative0 d0 = calculus.derivative<0, PRIMAL>();
    const Calculus::PrimalDerivative1 d1 = calculus.derivative<1, PRIMAL>();
    const Calculus::DualDerivative1 d1p = calculus.derivative<1, DUAL>();
    const Calculus::DualDerivative2 d2p = calculus.derivative<2, DUAL>();
    const Calculus::PrimalHodge1 h1 = calculus.primalHodge<1>();
    const Calculus::PrimalHodge2 h2 = calculus.primalHodge<2>();
    const Calculus::DualHodge2 h2p = calculus.dualHodge<2>();
    const Calculus::DualHodge3 h3p = calculus.dualHodge<3>();
    const LinearOperator<Calculus, 1, PRIMAL, 0, PRIMAL> ad1 = h3p * d2p * h1;
    const LinearOperator<Calculus, 2, PRIMAL, 1, PRIMAL> ad2 = h2p * d1p * h2;

	typedef EigenLinearAlgebraBackend::SolverSparseQR LinearAlgebraSolver;
	typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 0, PRIMAL, 0, PRIMAL> Solver;
	const Calculus::PrimalForm0 input_one_form_anti_derivated = ad1 * input_one_form;
    const Calculus::PrimalForm2 input_one_form_derivated = d1 * input_one_form;

	Display3DFactory<>::draw(viewer, calculus.sharp(input_one_form));
	Solver solver;
	solver.compute(ad1 * d0);
	Calculus::PrimalForm0  solution_curl_free = solver.solve(input_one_form_anti_derivated);

	typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 2, PRIMAL,2, PRIMAL> Solver2;
	Solver2 solver2;
	solver2.compute(d1 * ad2);
	Calculus::PrimalForm2 solution_div_free = solver2.solve(input_one_form_derivated);

	Calculus::PrimalForm1 solution_harmonic = input_one_form - d0*solution_curl_free - ad2*solution_div_free;
	Calculus::PrimalVectorField vfield = calculus.sharp(solution_harmonic);
	Display3DFactory<>::draw(viewer, vfield);*/
}

///////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv )
{
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (skeleton)")
		("output,o", po::value<std::string>(), "sliced vol file with orthogonal planes")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		("smallRadius,r", po::value<int>()->default_value(3), "small radius")
		("bigRadius,R", po::value<int>()->default_value(5), "big radius")
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
	double R = vm["bigRadius"].as<int>();
	double r = vm["smallRadius"].as<int>();
	
	

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

	trace.info() << "File             = " << inputFilename << std::endl;
	trace.info() << "Min image thres. = " << thresholdMin << std::endl;
	trace.info() << "Max image thres. = " << thresholdMax << std::endl;

	trace.info() << "Big radius     R = " << R << std::endl;

	trace.info() << "Small radius   r = " << r << std::endl;
	const double trivial_r = 3;
	trace.info() << "Trivial radius t = " << trivial_r << std::endl; // for orienting the directions given by the tensor.
	const double T = 0.1;
	trace.info() << "Feature thres. T = " << T << std::endl; // threshold for displaying features as red.

	const double size = 1.0; // size of displayed normals.

	KSpace ks;
	// Reads the volume
	trace.beginBlock( "Loading image into memory and build digital surface." );
	Image image = GenericReader<Image>::import(inputFilename );
	Z3i::DigitalSet aSet(image.domain());


	  ThresholdedImage thresholdedImage( image, thresholdMin, thresholdMax );
	trace.endBlock();
	trace.beginBlock( "Extracting boundary by scanning the space. " );
	ks.init( image.domain().lowerBound(),
			 image.domain().upperBound(), true );
	SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
	Surfel bel = Surfaces<KSpace>::findABel( ks, thresholdedImage, 100000 );
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
							  chi, trivial_r, l2, true );
	//! [DVCM3D-instantiation]

	trace.beginBlock( "Displaying VCM" );
	QApplication application(argc,argv);
	Viewer3D<> viewer( ks );
	Cell dummy;
	viewer.setWindowTitle("3D VCM viewer");
	//viewer << SetMode3D( dummy.className(), "Illustration" );
	viewer.show();
	SetFromImage<Z3i::DigitalSet>::append(aSet, image, thresholdMin, thresholdMax);
	
	DEC(aSet, image.domain(),viewer);

	/*typedef EigenDecomposition<3,double> LinearAlgebraTool;
	typedef LinearAlgebraTool::Matrix Matrix;
	vcm_surface.radiusTrivial();
	GradientColorMap<double> grad( 0, T );
	grad.addColor( Color( 128, 128, 255 ) );
	grad.addColor( Color( 255, 255, 255 ) );
	grad.addColor( Color( 255, 255, 0 ) );
	grad.addColor( Color( 255, 0, 0 ) );
	RealVector lambda; // eigenvalues of chi-vcm
	Matrix m;
	vector<Point> points;
	int i = 0;
	for ( S2NConstIterator it = vcm_surface.mapSurfel2Normals().begin(), 
			  itE = vcm_surface.mapSurfel2Normals().end(); it != itE; ++it )
    {
		Surfel s = it->first;
 
		Point kp = ks.sKCoords( s );
		RealPoint rp( 0.5 * (double) kp[ 0 ], 0.5 * (double) kp[ 1 ], 0.5 * (double) kp[ 2 ] );
		vcm_surface.getChiVCMEigenStructure( lambda, m, s );
		RealVector n = m.column(1);
		double ratio = lambda[ 1 ] / ( lambda[ 0 ] + lambda[ 1 ] + lambda[ 2 ] ); 
		if (i == 3000) {
			RealPoint theOtherPoint = rp;
			map<Point, VCMOnSurface::EigenStructure> mapPE = vcm_surface.mapPoint2ChiVCM();
			for (int c = 0; c < 100; c++) {
				cout << theOtherPoint << endl;
				Matrix m = mapPE[theOtherPoint].vectors;
				RealVector savedn = n;
				n = m.column(1);				
				cout << n << endl;
				theOtherPoint += n;
				points.push_back(theOtherPoint);
			}
		}
		i++;
		viewer.setFillColor( grad( ratio > T ? T : ratio ) );
	  	viewer << ks.unsigns( s );
		n *= size;
		viewer.setLineColor( Color::Black );
		viewer.addLine( rp + n, rp - n, 0.1 );
	}
	cout << points.size() << endl;
	for (const Point & p : points) {
		viewer << CustomColors3D(Color::Blue, Color::Blue) << p;
		}*/
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
