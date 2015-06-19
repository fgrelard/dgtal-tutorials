
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
#include <QtGui/qapplication.h>

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
#include "DGtal/geometry/surfaces/ChordGenericStandardPlaneComputer.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/shapes/parametric/Ball3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/writers/VolWriter.h"
#include <miniball/Seb.h>
#include <miniball/Seb_debug.h>
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/geometry/surfaces/estimation/VoronoiCovarianceMeasureOnDigitalSurface.h"

#include "geometry/DistanceToPointFunctor.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


template <typename Point>
class WeightedPoint {
public:
	WeightedPoint( const Point& _p, float _d ) : p(_p), d(_d) {}
	friend bool operator<(const WeightedPoint& it, const WeightedPoint& other) {
		if (it.d >= other.d) {
			return true;
		}
		return false;
	}
public:
	Point p;
	float d;
};


template <typename ImageFct, typename WeightedPoint, typename Point>
void checkPointForMedialAxis(ImageFct& imageFct, vector<WeightedPoint>& vPoints , const Point& p, float d) {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space,2> Distance;
	typedef DistanceToPointFunctor<Distance>         DistanceToPoint;
	typedef MetricAdjacency<Z3i::Space, 3>                Graph;
	typedef DistanceBreadthFirstVisitor< Graph, DistanceToPoint, std::set<Point> >
		DistanceVisitor;
	typedef typename DistanceVisitor::Node MyNode;
	typedef typename DistanceVisitor::Scalar MySize;	
	Graph             graph;
	DistanceToPoint   d2pfct( Distance(), p );
	DistanceVisitor   visitor( graph, d2pfct, p );
//trace.info() << "last= " << last << endl;
	MyNode node;
	bool add = true;
	visitor.expand(); //Go to the first neighbour	
	while ( ! visitor.finished() )
	{
	
		node = visitor.current(); // all the vertices of the same layer have been processed.
		
		if (node.second > 1) break; // we want to analyse only the direct neighbourhood (4- or 6- connexity)
		if (imageFct.domain().isInside(node.first) && node.second == 1) { //is inside domain
			float distanceToBoundary = imageFct(node.first);
			float minEnclosingRadius = sqrt(1+pow(d, 2));
			if (d <= 1 || distanceToBoundary >= minEnclosingRadius) {
				add = false;
			}
			visitor.expand();
		}
		else
			visitor.ignore();
	}
	if (add) {
		vPoints.push_back(WeightedPoint(p, d));
	}
}


template <typename Point, typename WeightedPoint>
void checkPointForMedialAxisNaive( vector<WeightedPoint>& vPoints, const Point& p, float d )
{
	typedef Z3i::Space Space;
	typedef Ball3D<Space> Ball3d;
	
	Ball3d ball(p[0], p[1], p[2], d);
	Point currentBallLowerBound = ball.getLowerBound();
	Point currentBallUpperBound = ball.getUpperBound();
	bool add = true;
	for (typename vector<WeightedPoint>::iterator it = vPoints.begin(), itE = vPoints.end(); it != itE; ++it) {
		Ball3d otherBall((*it).p[0], (*it).p[1], (*it).p[2], (*it).d);
		Point otherBallLowerBound = otherBall.getLowerBound();
		Point otherBallUpperBound = otherBall.getUpperBound();
		if (currentBallLowerBound >= otherBallLowerBound &&
			currentBallUpperBound <= otherBallUpperBound) {
			add = false;
		}
	}
	if (add) {
		vPoints.push_back(WeightedPoint( p, d ));
	}
}


template <typename Point, typename DGtalPoint>
void smallestEnclosingBall(set<DGtalPoint>& points, vector<Point>& sebPoints) {
	
	for (typename set<DGtalPoint>::iterator it = points.begin(),
			 itE = points.end();
		 it != itE;
		 ++it) {
		Point p(DGtalPoint::dimension, (*it).begin());
		sebPoints.push_back(p);
	}
}

template <typename ImageFct, typename WeightedPoint, typename Point>
void checkPointForLambdaMedialAxis(const ImageFct& imageFct, vector<WeightedPoint>& vPoints, const Point& p, float d, double lambda)  {
	typedef ExactPredicateLpSeparableMetric<Z3i::Space,2> Distance;
	typedef DistanceToPointFunctor<Distance>         DistanceToPoint;
	typedef MetricAdjacency<Z3i::Space, 1>                Graph;
	typedef DistanceBreadthFirstVisitor< Graph, DistanceToPoint, std::set<Point> >
		DistanceVisitor;
	typedef typename DistanceVisitor::Node MyNode;
	typedef typename DistanceVisitor::Scalar MySize;

	//smallest enclosing ball
	typedef typename Seb::Smallest_enclosing_ball<double> Miniball;
	typedef Seb::Point<double> SebPoint;
	
	Graph             graph;
	DistanceToPoint   d2pfct( Distance(), p );
	DistanceVisitor   visitor( graph, d2pfct, p );
	//trace.info() << "last= " << last << endl;
	MyNode node;

	set<Point> ballPoints;
	visitor.expand(); //Go to the first neighbour
	
	while ( ! visitor.finished() )
	{
		node = visitor.current(); // all the vertices of the same layer have been processed.
		
		if (node.second > 1) break; // we want to analyse only the direct neighbourhood (4- or 6- connexity) 
		if (imageFct.domain().isInside(node.first) && node.second == 1) { //is inside domain
			float distanceToBoundary = imageFct(node.first);
			if ( d > 1 && distanceToBoundary < d ) { 
				//Extended projection
				Point q = imageFct.getVoronoiVector(node.first);
				ballPoints.insert(q);
			}
 			visitor.expand();
		}
		else {
			visitor.ignore();
		}
	}
	vector<SebPoint> sebPoints;
	smallestEnclosingBall(ballPoints, sebPoints);
	if (sebPoints.size() > 0) {
		Miniball ball(Point::dimension, sebPoints);
		if (ball.radius() >= lambda * d) {
			vPoints.push_back(WeightedPoint(p, d));
		}
	}
   
	
}

template <typename Point>
void computeVCM(Viewer3D<> &viewer, vector<Point> & points, float R, float r, int T) {
	typedef typename vector<Point>::iterator Iterator;
	typedef Z3i::Space Space;
	typedef HyperRectDomain<Space> Domain;
	typedef ExactPredicateLpSeparableMetric<Space, 2> Metric;          // L2-metric type
	typedef functors::HatPointFunction<Point,double>  KernelFunction;  // chi function type 
	typedef EigenDecomposition<3,double> LinearAlgebraTool;
	typedef LinearAlgebraTool::Matrix Matrix;
	typedef Z3i::KSpace KSpace;
	typedef VoronoiCovarianceMeasure<Space,Metric> VCM;
	typedef functors::HatPointFunction<Point,double> KernelFunction;
	typedef Z3i::RealPoint RealPoint;
	typedef Z3i::RealVector RealVector;
	typedef KSpace::Cell Cell;
	typedef ChordNaivePlaneComputer<Z3i::Space, Z3i::Point, int32_t> PlaneComputer;

	Iterator it = points.begin();
    Iterator itE = points.end();
	
	Metric l2;
	VCM vcm(R,r, l2, true);
	vcm.init( points.begin(), points.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, r );
	//ks.init( image.domain().lowerBound(),
	//		 image.domain().upperBound(), true );
	// Flat zones are metallic blue, slightly curved zones are white,
	// more curved zones are yellow till red.
	GradientColorMap<double> colormap( 0.0, T );
	colormap.addColor( Color( 128, 128, 255,255 ) );
	colormap.addColor( Color( 255, 255, 255,255 ) );
	colormap.addColor( Color( 255, 255, 0, 255 ) );
	colormap.addColor( Color( 255, 0, 0, 255 ) );
	Matrix vcm_r, evec;
	RealVector eval;
	Cell dummy;
	double size = 1.0;
	for (; it != itE; ++it )
	{
		// Compute VCM and diagonalize it.
		vcm_r = vcm.measure( chi, *it );
		LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
		double feature = eval[ 1 ] / ( eval[ 0 ] +  eval[ 1 ] + eval[2] );
		// Display normal
		RealVector n = evec.column( 1 );
		RealVector n2 = evec.column( 0 );

		
		RealPoint rp( (*it)[ 0 ], (*it)[ 1 ], (*it)[2] );
		Point p1((*it)[0] + n[0] * 2 + n2[0] * 2, (*it)[1] + n[1] * 2 + n2[1] * 2 , (*it)[2] + n[2] * 2 + n2[2] * 2);
		Point p2((*it)[0] + n[0] * -2 +  n2[0] * 2, (*it)[1] + n[1] * -2 + n2[1] * 2, (*it)[2] + n[2] * -2 + n2[2] * 2 );
		Point p3((*it)[0] + n[0] * 2 + n2[0] * -2 , (*it)[1] + n[1] * 2 + n2[1] * -2 , (*it)[2] + n[2] * 2 + n2[2] * -2 );
		Point p4((*it)[0] + n[0] * -2 + n2[0] * -2 , (*it)[1] + n[1] * -2 + n2[1] * -2 , (*it)[2] + n[2] * -2 + n2[2] * -2 );						
		//viewer.setFillColor( colormap( 0.0 ) );
		viewer.setFillColor(Color::Red);
		//viewer.addQuad(p1, p3, p4, p2);
		viewer.setFillColor( Color( 100, 100, 140, 150 ) );
		viewer << (*it);
		n *= size;
		viewer.setLineColor( Color::Black );
		viewer.addLine( rp + n, rp - n, 0.3 );
	}
}



int main( int argc, char** argv )
{
	using namespace DGtal;
	using namespace DGtal::Z3i;
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef ImageSelector<Z3i::Domain, float>::Type FloatImage;
	typedef functors::IntervalForegroundPredicate<Image> Binarizer;
	typedef WeightedPoint<Point> WeightedPoint;
	typedef functors::NotPointPredicate<Binarizer> NotPredicate;
	typedef  DistanceTransformation<Z3i::Space, Binarizer, Z3i::L2Metric> DTL2;


	if ( argc <= 2 ) return 1;

	QApplication app(argc, argv);
	Viewer3D<> viewer;
	Cell dummy;
	viewer << SetMode3D( dummy.className(), "Basic" );
	viewer.show();

	Image img  = VolReader<Image>::importVol( argv[ 1 ] );

	
	double lambda = atof(argv[2]);
	Z3i::Domain domain = img.domain();

	Binarizer binarizer(img, 0, 255);
	NotPredicate nBinarizer(binarizer);
	
	trace.beginBlock("Computing distance map");
    DTL2 dtL2(&domain, &binarizer, &Z3i::l2Metric);
	trace.endBlock();

	trace.beginBlock("Computing medial axis");
	set<WeightedPoint> pSet;
	for (DTL2::Domain::Iterator it = domain.begin(), itE = domain.end(); it != itE; ++it) {
		if (img(*it) >= 1)
			pSet.insert( WeightedPoint( (*it), dtL2(*it) ) );
	}
	vector<WeightedPoint> vPoints;
	vPoints.push_back( *(pSet.begin()) );
	set<WeightedPoint>::iterator it = pSet.begin();
	++it;

	int nb = pSet.size();
	int i = 0;
	for (set<WeightedPoint>::iterator itE = pSet.end(); it != itE; ++it, ++i) {
	    checkPointForMedialAxis(dtL2, vPoints, (*it).p, (*it).d);
		trace.progressBar(i, nb);
	}
	trace.endBlock();

	trace.beginBlock("Displaying");

	Board3D<Space, KSpace> board;
	//  board << SetMode3D( d2.domain().className(), "Paving" );
	
	board << SetMode3D(domain.className(), "Paving");
	Color color(100,100,140,150);
	vector<Point> points;
	for ( vector<WeightedPoint>::iterator it = vPoints.begin(), itE = vPoints.end(); it != itE; ++it) 
	{
	  	viewer  << CustomColors3D(color, color) << (*it).p;
		points.push_back( (*it).p );
	}
//	computeVCM(viewer, points, 5, 3, 0.3);

	Image outImage(img.domain());
	DGtal::imageFromRangeAndValue(points.begin(), points.end(), outImage, 10);
	VolWriter<Image>::exportVol("/home/florent/test_img/medialaxis.vol", outImage);
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
	std::cout << endl;
	trace.endBlock();
	return 0;
}

