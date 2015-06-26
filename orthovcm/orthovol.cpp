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
#include "DGtal/geometry/volumes/estimation/VoronoiCovarianceMeasure.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

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
  typedef ImageContainerBySTLVector<Domain, unsigned char> Image;
  typedef VoronoiCovarianceMeasure<Space,Metric> VCM;
  typedef functors::HatPointFunction<Point,double> KernelFunction;


  string inputSDP = argv[1];
  std::vector<Point> pts = PointListReader<Point>::getPointsFromFile(inputSDP); 
  const double R = atof(argv[2]);
  trace.info() << "Big radius   R = " << R << std::endl;
  const double r = atof(argv[3]);
  trace.info() << "Small radius r = " << r << std::endl;
  //  const double T = 0.1;
  //  trace.info() << "Feature thres. T = " << T << std::endl; // threshold for displaying features as red.
  const double size = 35.0; // size of displayed normals
  string inputVOL = argv[4];
  string outputFName = argv[5];
  string sliceVOL = argv[6];
  Image image = VolReader<Image>::importVol(inputVOL);
  Image sliceImageVol = VolReader<Image>::importVol(sliceVOL);
  std::vector<Point> ptsvol;
  for(Domain::ConstIterator it=image.domain().begin(); it != image.domain().end(); ++it){
    if(image(*it) >= 1){
      ptsvol.push_back(*it);
     
    }
  }
  
  Metric l2;
  VCM vcm( R, ceil( r ), l2, true );
  vcm.init( ptsvol.begin(), ptsvol.end() );
  Domain domain = vcm.domain();
  KernelFunction chi( 1.0, r );

  // Flat zones are metallic blue, slightly curved zones are white,
  // more curved zones are yellow till red.
  Matrix vcm_r, evec;
  RealVector eval;
  Viewer3D<> viewer;
  viewer.show();

  typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;
  const int IMAGE_PATCH_WIDTH = 100;
  unsigned int sliceNumber = 0;
  Z3i::Domain domain3Dyup(image.domain().lowerBound() + Z3i::Point(-IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH), image.domain().upperBound() + Z3i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
  DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
  DGtal::functors::Identity idV;

  // for ( std::vector<Point>::const_iterator it = ptsvol.begin(), itE = ptsvol.end();
  //       it != itE; ++it )
  //   {
  //     // Compute VCM and diagonalize it.
  //     viewer.setFillColor(Color::Red);
  //     viewer.setFillTransparency(25);

  //     viewer << *it;
  //   }
  int i = 0;
  	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
  for ( std::vector<Point>::const_iterator it = pts.begin(), itE = pts.end();
        it != itE; ++it )
    {
      // Compute VCM and diagonalize it.
		viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<<*it;
      vcm_r = vcm.measure( chi, *it );
      LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
      //double feature = eval[1] /  (eval[0] + eval[ 1 ]+  eval[3]);
      // // Display normal
      RealVector n = evec.column( 2 );
      n*=size;
      RealPoint p( (*it)[ 0 ], (*it)[ 1 ], (*it)[2] );
      RealVector n2 = evec.column( 1 );
      n2*=size;
      //std::cout << evec.column(0) << std::endl;      
      DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, p, evec.column(0), IMAGE_PATCH_WIDTH);
      ImageAdapterExtractor extractedImage(sliceImageVol, domainImage2D, embedder, idV);
      std::string outName;
      outName += outputFName + "_" + std::to_string(sliceNumber) + ".pgm";
      sliceNumber++;
      GenericWriter<ImageAdapterExtractor>::exportFile(outName, extractedImage);

      viewer.setLineColor(Color::Red);
	  viewer.addLine(p+n, p-n, 0.1);
	  viewer.setLineColor(Color::Green);
	  viewer.addLine(p+n2, p-n2, 0.1);
      viewer.setLineColor(Color::Red);
      viewer.setFillColor(Color::Red);
      viewer.setFillTransparency(150);
	  if (i%20 == 0) {
		  //viewer.addQuad(p-n-n2,p-n+n2,p+n+n2,p+n-n2);
	  }
	  i++;
    }

  for (auto it = sliceImageVol.domain().begin(), itE = sliceImageVol.domain().end(); it!=itE; ++it) {
	  if (sliceImageVol(*it) > 0) {
		  viewer << CustomColors3D(Color::Blue, Color::Blue) << *it;
	  }
  }
  viewer << Viewer3D<>::updateDisplay;
  application.exec();
  return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
