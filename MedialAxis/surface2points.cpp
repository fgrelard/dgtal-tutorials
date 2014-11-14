
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
 * @file surface2points.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/05/02
 *
 * Computes the points of some digital surface.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/images/ImageSelector.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


template <typename Distance>
struct DistanceToPointFunctor {

  typedef typename Distance::Space Space;
  typedef typename Distance::Value Value;
  typedef typename Space::Point    Point;

  Point p;
  DistanceToPointFunctor( Clone<Distance> distance,
                          const Point& aP )
    : myDistance( distance ), p( aP ) {}

  Value operator()( const Point& q ) const
  {
    return myDistance( p, q );
  }
  Distance myDistance;
};

// A measure is a function
template <typename ImageFct>
class DistanceToMeasure {
public:
  typedef typename ImageFct::Value  Value;
  typedef typename ImageFct::Point  Point;
  typedef typename ImageFct::Domain Domain;
  typedef typename Domain::Space    Space;

public:
  DistanceToMeasure( const ImageFct& measure, Value rmax = 10.0 )
    : myMeasure( measure ), myDistance2( myMeasure.domain() ),
      myR2Max( rmax*rmax )
  {
    init( myMeasure );
  }

  void init( const ImageFct& measure )
  {
    double       nb = myDistance2.domain().size();
    unsigned int i  = 0;
    trace.progressBar( i, nb );

    for ( typename Domain::ConstIterator it = myDistance2.domain().begin(),
            itE = myDistance2.domain().end(); it != itE; ++it, ++i )
      {
        if ( ( i % 100 ) == 0 ) trace.progressBar( i, nb );
        myDistance2.setValue( *it, computeDistance2( *it ) );
      }
  }

  Value distance2( const Point& p ) const
  {
    return myDistance2( p );
  }

  Value computeDistance2( const Point& p )
  {
    typedef ExactPredicateLpSeparableMetric<Space,2> Distance;
    typedef DistanceToPointFunctor<Distance>         DistanceToPoint;
    typedef MetricAdjacency<Space, 1>                Graph;
    typedef DistanceBreadthFirstVisitor< Graph, DistanceToPoint, std::set<Point> >
      DistanceVisitor;
    typedef typename DistanceVisitor::Node MyNode;
    typedef typename DistanceVisitor::Scalar MySize;

    Value             m  = NumberTraits<Value>::ZERO;
    Value             d2 = NumberTraits<Value>::ZERO;
    Graph             graph;
    DistanceToPoint   d2pfct( Distance(), p );
    DistanceVisitor   visitor( graph, d2pfct, p );


    Value last = d2pfct( p );
    MyNode node;
    while ( ! visitor.finished() )
    {
      node = visitor.current();
      if ( ( node.second != last )) break; // all the vertices of the same layer have been processed.
      if ( node.second > myR2Max ) { d2 = m * myR2Max; break; }
      if ( myMeasure.domain().isInside( node.first ) )
	  {
		  last = node.second;
          visitor.expand();
	  }
      else
        visitor.ignore();
    }
    return d2 / m;
  }

public:
  const ImageFct& myMeasure;
  ImageFct myDistance2;
  Value myR2Max;
};

int main( int argc, char** argv )
{
  using namespace DGtal;
  using namespace DGtal::Z3i;
  typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
  typedef ImageSelector<Z3i::Domain, float>::Type FloatImage;
  typedef DistanceToMeasure<FloatImage>                 Distance;
  typedef functors::IntervalForegroundPredicate<Image> Binarizer;
  typedef  DistanceTransformation<Z3i::Space, Binarizer, Z3i::L2Metric> DTL2;
  if ( argc <= 2 ) return 1;
  Image img  = VolReader<Image>::importVol( argv[ 1 ] );
  Z3i::Domain domain = img.domain();
  double           rmax = atof( argv[ 2 ] );


  Binarizer binarizer(img, 0, 135);
  
  trace.beginBlock("Computing distance map");
  DTL2 dtL2(&domain, &binarizer, &Z3i::l2Metric);
  trace.endBlock();
  
  FloatImage     fimg( dtL2.domain() );
  FloatImage::Iterator outIt = fimg.begin();
  unsigned int max = 0;
  for(DTL2::ConstRange::ConstIterator it = dtL2.constRange().begin(), 
			itend=dtL2.constRange().end();
		it!=itend;
		++it)
    {
		if( (*it) > max ) 
			max=(*it);
    }

  
  
  for ( Image::ConstIterator it = img.begin(), itE = img.end();
        it != itE; ++it )
    {
      float v = ((float)*it) / max;
      *outIt++ = v;
    }
  trace.beginBlock( "Computing delta-distance." );
  Distance     delta( fimg, rmax );
  const FloatImage& d2 = delta.myDistance2;
  trace.endBlock();

  float m = 0.0;
  for ( Domain::ConstIterator it = d2.domain().begin(),
          itE = d2.domain().end(); it != itE; ++it )
    {
      Point p = *it;
      float v = sqrt( d2( p ) );
      m = std::max( v, m );
    }

  GradientColorMap<float> cmap_grad( 0,  m);
  cmap_grad.addColor( Color( 255, 255, 255 ) );
  cmap_grad.addColor( Color( 255, 255, 0 ) );
  cmap_grad.addColor( Color( 255, 0, 0 ) );
  cmap_grad.addColor( Color( 0, 255, 0 ) );
  cmap_grad.addColor( Color( 0,   0, 255 ) );
  cmap_grad.addColor( Color( 0,   0, 0 ) );
  Board3D<Space, KSpace> board;
  //  board << SetMode3D( d2.domain().className(), "Paving" );


  for ( typename Domain::ConstIterator it = d2.domain().begin(),
          itE = d2.domain().end(); it != itE; ++it )
    {
      Point p = *it;
      float v = sqrt( d2( p ) );
      v = std::max( v, 0.0f );
      board  << p;
    }
  std::cout << endl;
  board.saveOBJ("delta2.obj");
  return 0;
}

