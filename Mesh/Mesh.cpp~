#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "VoronoiCovarianceMeasure.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/boards/Board2D.h"


///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////


int main( int /* argc */, char** /* argv */ )
{
  const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
  typedef Z2i::Space Space;
  typedef Z2i::Point Point;
  typedef Z2i::RealPoint RealPoint;
  typedef Z2i::RealVector RealVector;
  typedef HyperRectDomain<Space> Domain;
  typedef DGtal::Z2i::L2Metric Metric; // L2-metric
  typedef EigenDecomposition<2,double> LinearAlgebraTool;
  typedef LinearAlgebraTool::Matrix Matrix;
  
  typedef VoronoiCovarianceMeasure<Z2i::Space,Metric> VCM;
  typedef functors::HatPointFunction<Point,double> KernelFunction;
  // Gets the points
  vector<unsigned int> vPos;
  vPos.push_back(0);
  vPos.push_back(1);
  string inputSDP = examplesPath + "flower-30-8-3.sdp";
  // string inputSDP = examplesPath + "samples/ellipse-20-7-0.4.sdp";
  // string inputSDP = examplesPath + "samples/accflower-20-5-5-0.1.sdp";
  trace.info() << "Reading input 2d discrete points file: " << inputSDP; 
  std::vector<Point> pts = PointListReader<Point>::getPointsFromFile(inputSDP, vPos); 
  trace.info() << " [done] " << std::endl ; 
  const double R = 20;
  trace.info() << "Big radius   R = " << R << std::endl;
  const double r = 5;
  trace.info() << "Small radius r = " << r << std::endl;
  const double T = 0.1;
  trace.info() << "Feature thres. T = " << T << std::endl; // threshold for displaying features as red.
  const double size = 3.0; // size of displayed normals

  Metric l2;
  VCM vcm( 20.0, 5.0, l2, true );
  trace.info() << "VCM created " << endl;
  vcm.init( pts.begin(), pts.end() );
  trace.info() << "VCM computed " << endl;
  Domain domain = vcm.domain();
  KernelFunction chi( 1.0, r );
  // Flat zones are metallic blue, slightly curved zones are white,
  // more curved zones are yellow till red.
  GradientColorMap<double> colormap( 0.0, T );
  colormap.addColor( Color( 128, 128, 255 ) );
  colormap.addColor( Color( 255, 255, 255 ) );
  colormap.addColor( Color( 255, 255, 0 ) );
  colormap.addColor( Color( 255, 0, 0 ) );
  Board2D board;
  Matrix vcm_r, evec;
  RealVector eval;
  for ( std::vector<Point>::const_iterator it = pts.begin(), itE = pts.end();
        it != itE; ++it )
    {
      // Compute VCM and diagonalize it.
      vcm_r = vcm.measure( chi, *it );
      LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
      double feature = eval[ 0 ] / ( eval[ 0 ] +  eval[ 1 ] );
      board << CustomStyle( it->className(), 
                            new CustomColors( Color::Black,  colormap( feature > T ? T : feature ) ) )
            << *it;
      // Display normal
      RealVector normal = evec.column( 1 );
      RealPoint p( (*it)[ 0 ], (*it)[ 1 ] ); 
      Display2DFactory::draw( board, size*normal, p );
      Display2DFactory::draw( board, -size*normal, p );
    }      
  board.saveSVG("dvcm-hat-r.svg");
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
