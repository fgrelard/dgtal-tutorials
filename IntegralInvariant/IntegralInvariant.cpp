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
 * @file exampleIntegralInvariantCurvature3D.cpp
 * @ingroup Examples
 * @author Jérémy Levallois (\c jeremy.levallois@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), INSA-Lyon, France
 * LAboratoire de MAthématiques - LAMA (CNRS, UMR 5127), Université de Savoie, France
 *
 * @date 2012/12/17
 *
 * An example file named exampleIntegralInvariantCurvature3D.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"

// Shape construction
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/SimpleThresholdForegroundPredicate.h"
#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
/// Estimator
#include "DGtal/geometry/surfaces/estimation/IIGeometricFunctors.h"
#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"
#include "DGtal/kernel/BasicPointPredicates.h"
// Drawing
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include <QtGui/QApplication>
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/base/CountedConstPtrOrConstPtr.h"
#include "IntegralInvariantVolumeEstimatorAdapted.h"

using namespace DGtal;
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

template <typename T, typename S>
class Couple {
public:
	Couple(const T& f, const T& s, float _d) : first(f), end(s), d(_d) {}
public:
	T first;
	T end;
	float d;
};

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
	typedef MetricAdjacency<Z3i::Space, 1>                Graph;
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

template <typename WeightedSurfel, typename SurfelConstIterator, typename Range, typename WeightedPoint,
		  typename KSpace>
void closestMedialAxisPoint(vector<WeightedSurfel> & vIterators, const Range & range,
							const vector<WeightedPoint> & vPoints, KSpace KSpaceShape) {
	for (SurfelConstIterator it = range.begin(); it != range.end(); ++it) {
		double maximum = INT_MAX;
		Z3i::Point surfPoint = KSpaceShape.sCoords(*it);
		for (typename vector<WeightedPoint>::const_iterator it2 = vPoints.begin(); it2 != vPoints.end(); ++it2) {
			double euclideanDistance = sqrt(pow((surfPoint[0] - (*it2).p[0]), 2) +
											pow((surfPoint[1] - (*it2).p[1]), 2) +
											pow((surfPoint[2] - (*it2).p[2]), 2));
			if (euclideanDistance < maximum) {
				maximum = euclideanDistance;
			}
		}
		WeightedSurfel w(*it, maximum);
		vIterators.push_back(w);
	}
}

///////////////////////////////////////////////////////////////////////////////

using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
    if ( argc != 5 )
    {
        trace.error() << "Usage: " << argv[0]
                               << " <fileName.vol> <threshold> <radius>" << std::endl;
        trace.error() << "Example : "<< argv[0] << " Al.150.vol 0 7" << std::endl;
        return 0;
    }

    trace.beginBlock ( "Example IntegralInvariantCurvature3D" );
    trace.info() << "Args:";
    for ( int i = 0; i < argc; ++i )
        trace.info() << " " << argv[ i ];
    trace.info() << endl;

    double h = 1.0;
    unsigned int threshold = std::atoi( argv[ 2 ] );
	double leeway = std::atof(argv[3]);
	double scaling_factor = std::atof(argv[4]);
    /// Construction of the shape from vol file
    typedef ImageSelector< Z3i::Domain, bool >::Type Image;
    typedef functors::SimpleThresholdForegroundPredicate< Image > ImagePredicate;
    typedef LightImplicitDigitalSurface< Z3i::KSpace, ImagePredicate > MyLightImplicitDigitalSurface;
    typedef DigitalSurface< MyLightImplicitDigitalSurface > MyDigitalSurface;
	typedef ImageSelector<Z3i::Domain, float>::Type FloatImage;
	typedef functors::IntervalForegroundPredicate<Image> Binarizer;
	typedef functors::NotPointPredicate<Binarizer> NotPredicate;
	typedef DistanceTransformation<Z3i::Space, Binarizer, Z3i::L2Metric> DTL2;
	typedef DepthFirstVisitor< MyDigitalSurface > Visitor;
    typedef GraphVisitorRange< Visitor > VisitorRange;
    typedef VisitorRange::ConstIterator SurfelConstIterator;
	typedef VisitorRange::NodeConstIterator NodeConstIterator;
	typedef VisitorRange::Vertex Surfel;
	typedef WeightedPoint<Surfel> WeightedSurfel;
	typedef WeightedPoint<Z3i::Point> WeightedPoint;

	typedef Couple<SurfelConstIterator, SurfelConstIterator> couple;
	
    std::string filename = argv[1];
    Image image = VolReader<Image>::importVol( filename );
    ImagePredicate predicate = ImagePredicate( image, threshold );

    Z3i::Domain domain = image.domain();

    Z3i::KSpace KSpaceShape;

    bool space_ok = KSpaceShape.init( domain.lowerBound(), domain.upperBound(), true );
    if (!space_ok)
    {
      trace.error() << "Error in the Khalimsky space construction."<<std::endl;
      return 2;
    }
	Binarizer binarizer(image, 0, 255);
	NotPredicate nBinarizer(binarizer);
	
	trace.beginBlock("Computing distance map");
    DTL2 dtL2(&domain, &binarizer, &Z3i::l2Metric);
	trace.endBlock();
	trace.beginBlock("Computing medial axis");
	set<WeightedPoint> pSet;
	for (DTL2::Domain::Iterator it = domain.begin(), itE = domain.end(); it != itE; ++it) {
		if (image(*it) >= 1)
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
		
	SurfelAdjacency< Z3i::KSpace::dimension > SAdj( true );
    Z3i::KSpace::Surfel bel = Surfaces< Z3i::KSpace >::findABel( KSpaceShape, predicate, 100000 );
    MyLightImplicitDigitalSurface LightImplDigSurf( KSpaceShape, predicate, SAdj, bel );
    MyDigitalSurface digSurf( LightImplDigSurf );
	VisitorRange range( new Visitor( digSurf, *digSurf.begin() ) );
   	vector<WeightedSurfel> vIterators;
		

	closestMedialAxisPoint<WeightedSurfel, SurfelConstIterator>(vIterators, range, vPoints, KSpaceShape);

	
    //! Integral Invariant stuff
    //! [IntegralInvariantUsage]


    typedef functors::IIMeanCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
    typedef IntegralInvariantVolumeEstimatorAdapted< Z3i::KSpace, ImagePredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
    
    // For computing Gaussian curvature instead, for example, change the two typedef above by : 
    // typedef functors::IIGaussianCurvature3DFunctor<Z3i::Space> MyIICurvatureFunctor;
    // typedef IntegralInvariantCovarianceEstimator< Z3i::KSpace, ImagePredicate, MyIICurvatureFunctor > MyIICurvatureEstimator;
    // and it's done. The following part is exactly the same.

    typedef MyIICurvatureFunctor::Value Value;
	std::vector< Value > results;
    std::back_insert_iterator< std::vector< Value > > resultsIt( results ); // output iterator for results of Integral Invariant curvature computation
	VisitorRange range2 ( new Visitor( digSurf, *digSurf.begin() ) );
	SurfelConstIterator abegin = range2.begin();
    SurfelConstIterator aend = range2.end();


	trace.beginBlock("Estimating curvature");
	MyIICurvatureFunctor curvatureFunctor; // Functor used to convert volume -> curvature
	float radius = 10;
	curvatureFunctor.init( h, radius ); // Initialisation for a grid step and a given Euclidean radius of convolution kernel
	MyIICurvatureEstimator curvatureEstimator( curvatureFunctor ); 
	curvatureEstimator.attach( KSpaceShape, predicate ); // Setting a KSpace and a predicate on the object to evaluate
	curvatureEstimator.setParams( radius / h ); // Setting the digital radius of the convolution kernel
	curvatureEstimator.init( h, abegin, aend ); // Initialisation for a given h, and a range of surfels
	curvatureEstimator.eval<Surfel, WeightedSurfel,  vector<Value>, SurfelConstIterator >( results, vIterators, leeway, scaling_factor ); // Computation
	trace.endBlock();
    
    //! [IntegralInvariantUsage]

    /// Drawing results
    Value min = numeric_limits < Value >::max();
    Value max = numeric_limits < Value >::min();
    for ( unsigned int i = 0; i < results.size(); ++i )
    {
        if ( results[ i ] < min )
        {
            min = results[ i ];
        }
        else if ( results[ i ] > max )
        {
            max = results[ i ];
        }
	}
	 
    QApplication application( argc, argv );
    typedef Viewer3D<Z3i::Space, Z3i::KSpace> Viewer;
    Viewer viewer( KSpaceShape );
    viewer.setWindowTitle("example Integral Invariant 3D");
    viewer.show();

    typedef GradientColorMap< Value > Gradient;
	trace.info() << min << " " << max;
	Gradient cmap_grad( min, max);
    cmap_grad.addColor( Color( 50, 50, 255 ) );
    cmap_grad.addColor( Color( 255, 0, 0 ) );
    cmap_grad.addColor( Color( 255, 255, 10 ) );

    VisitorRange range3( new Visitor( digSurf, *digSurf.begin() ) );
    abegin = range3.begin();

    Z3i::KSpace::Cell dummy_cell;
    viewer << SetMode3D( dummy_cell.className(), "Basic" );
	trace.info() << results.size() << " " << vIterators.size() << endl;
    for ( unsigned int i = 0; i < vIterators.size(); ++i )
    {
        viewer << CustomColors3D( Color::Black, cmap_grad( results[i] ))
               << KSpaceShape.unsigns( *abegin );
        ++abegin;
	}


	for (vector<WeightedPoint>::iterator it = vPoints.begin();
		 it != vPoints.end();
		 ++it) {
		viewer << CustomColors3D(Color::Black, Color(150,150,150,100)) << (*it).p;
	}
    viewer << Viewer3D<>::updateDisplay;

    trace.endBlock();
    return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
