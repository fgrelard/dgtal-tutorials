#include <iostream>
#include <cstdio>
#include <set>
#include <QtGui/qapplication.h>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/images/ImageSelector.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/shapes/parametric/Ball3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"

#include "DGtal/math/linalg/EigenSupport.h"
#include "DGtal/dec/DiscreteExteriorCalculus.h"
#include "DGtal/dec/DiscreteExteriorCalculusSolver.h"

#include "SurfaceUtils.h"
#include "SurfacePoint.h"
#include "../CreateCurve/Distance.h"
#include "../CreateCurve/Ball.h"
#include "PointUtil.h"
#include "../Tangent3D/SliceUtils.h"
#include "Path.h"

using namespace std;
using namespace DGtal;
using namespace Z3i;
namespace po = boost::program_options;


template <typename Vector>
bool allSigns(const Vector& v1, const Vector& v2) {
	double x = v1[0] * v2[0];
	double y = v1[1] * v2[1];
	double z = v1[2] * v2[2];
	return (x <= 0 && y <= 0 && z <= 0);
}


/**
 * Returns the point contained in surfacePointSet closest to pSkeleton
 */
template <typename Point>
Point pointCorrespondingToMinDistance(const Point& pSkeleton, const set<Point>& surfacePointSet) {
	Point minimumDistPoint = *(std::min_element(surfacePointSet.begin(), surfacePointSet.end(), [&](const Point& one, const Point& two) {
				return euclideanDistance(one, pSkeleton) < euclideanDistance(two, pSkeleton);
			}));
	return minimumDistPoint;
}

/**
 * Checks if curve is closed (all points have 2 neighbors)
 */
template <typename Point>
bool closedCurve(const set<Point>& points) {
	int cpt = 0;
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		int previouscpt = cpt;
		for (auto secondIt = points.begin(), secondIte = points.end(); secondIt != secondIte; ++secondIt) {
			if (PointUtil::areAlmostSimilar(*it, *secondIt) && *secondIt != *it) {
				cpt++;
			}
		}
		if (cpt - previouscpt < 2 || cpt - previouscpt > 3) cpt = 0;
	}
	return cpt >= points.size() * 2;
}



template <typename Point>
map<Point, Point> computeDSSOnPath(const vector<Point>& path) {
	typedef StandardDSS6Computer<typename vector<Point>::const_iterator,int,8> SegmentComputer;
	typedef GreedySegmentation<SegmentComputer> Segmentation;

	map<Point, Point> mapPointToDirection;
	SegmentComputer algo;
	Segmentation s(path.begin(), path.end(), algo);
	
	for (auto it = s.begin(), itE = s.end(); it!=itE; ++it) {
		SegmentComputer current(*it);
	    Point direction;
		RealPoint intercept, thickness;
		current.getParameters(direction, intercept, thickness);
		if (direction != Point::zero) {
			for (const Point& pointInDSS : current) {
				mapPointToDirection[pointInDSS] = direction;
			}
		}
	}
	return mapPointToDirection;
}


template <typename Point>
Point computeProjection(const Point& initialPoint, const Point& pivot, const DigitalSet& points) {
	Point toPivot = pivot - initialPoint;
	Point projected = pivot + toPivot;
	double mini = std::numeric_limits<double>::max();
	Point p;
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		double distanceToProjected = euclideanDistance(*it, projected);
		if (distanceToProjected < mini) {
			mini = distanceToProjected;
			p = *it;
		}
	}
	return p;
}

template <typename Space, typename KSpace, typename TCalculus, DGtal::Order order, DGtal::Duality duality>
void
draw(Display3D<Space, KSpace>& display, const DGtal::KForm<TCalculus, order, duality>& kform, double cmap_min, double cmap_max)
{
    BOOST_STATIC_ASSERT(( TCalculus::dimension == 3 ));
    ASSERT( kform.myCalculus );

    typedef typename TCalculus::Scalar Scalar;
    typedef typename TCalculus::SCell SCell;
    typedef typename TCalculus::Index Index;

    if (cmap_min == 0 && cmap_max == 0)
    {
        bool first = true;
        for (Index index=0; index<kform.myContainer.rows(); index++)
        {
            const Scalar value = kform.myContainer(index);
            if (!std::isfinite(abs(log10(value)))) continue;
            if (first || cmap_min > abs(log10(value))) cmap_min = abs(log10(value));
            if (first || cmap_max < abs(log10(value))) cmap_max = abs(log10(value));
            first = false;
        }
    }
	double denominator = cmap_max;
	cmap_min = 0.0;
	cmap_max = 1.0;
	
    typedef typename DGtal::HueShadeColorMap<Scalar, 2 > ColorMap;
    const ColorMap colormap(cmap_min, cmap_max);
                            
    for (typename TCalculus::Index index=0; index<kform.myCalculus->kFormLength(order, duality); index++)
    {
        const SCell& cell = kform.myCalculus->getSCell(order, duality, index);
        ASSERT(kform.myCalculus->myKSpace.sSign(cell) == TCalculus::KSpace::POS);

        const bool& flipped = kform.myCalculus->isSCellFlipped(cell);
        SCell displayed_cell = cell;
        if (flipped) displayed_cell = kform.myCalculus->myKSpace.sOpp(cell);

        Scalar displayed_value = kform.myContainer(index);
        if (flipped) displayed_value = -displayed_value;

        display << SetMode3D(cell.className(), "Basic");
        if (std::isfinite(displayed_value))
			display << DGtal::CustomColors3D(DGtal::Color::Black, colormap(abs(log10(displayed_value))/denominator) );
        else
            display << DGtal::CustomColors3D(DGtal::Color::Black, DGtal::Color::White);

        display << displayed_cell;
    }
}

template <typename Vertex, typename Result, typename Tracker, typename DirIterator>
map<Vertex, double> computeMapInOneDirection(Vertex& s,  const Result& result, Tracker& tracker, DirIterator& itDirs,double currentScalar ) {
	map<Vertex, double> aMapVertexToScalar;
	if ( tracker->adjacent( s, *itDirs, true ) ) {
		typename Result::Calculus::Index index = result.myCalculus->getSCellIndex(s);
		const bool& flipped = result.myCalculus->isSCellFlipped(s);
		typename Result::Calculus::DualForm1::Scalar otherScalar = result.myContainer(index);
		if (flipped) otherScalar = -otherScalar;
		if (currentScalar < otherScalar) aMapVertexToScalar[s] = otherScalar;
	}
	if ( tracker->adjacent( s, *itDirs, false ) ) {
		typename Result::Calculus::Index index = result.myCalculus->getSCellIndex(s);
		const bool& flipped = result.myCalculus->isSCellFlipped(s);
		typename Result::Calculus::DualForm1::Scalar otherScalar = result.myContainer(index);
		if (flipped) otherScalar = -otherScalar;
		if (currentScalar < otherScalar) aMapVertexToScalar[s] = otherScalar;
	}
	return aMapVertexToScalar;
}

template <typename Vertex, typename Boundary, typename Result>
bool belongsToBisector(const Vertex& vertex, const Result& result, const Boundary& boundary, double currentScalar) {
	auto itDirs =  result.myCalculus->myKSpace.sDirs(vertex);
	auto tracker = boundary.container().newTracker( vertex );
	Vertex s;
	for (; itDirs != 0; ++itDirs) {
	    map<Vertex, double> aMapVertexToScalar = computeMapInOneDirection(s, result, tracker, itDirs, currentScalar);
		if (aMapVertexToScalar.size() == 2) {
			// for (auto it = aMapVertexToScalar.begin(), ite = aMapVertexToScalar.end(); it != ite; ++it) {
			// 	auto trackerNeighbors = boundary.container().newTracker( it->first );
			// 	auto itOrthDirs =  result.myCalculus->myKSpace.sDirs( it->first );
			// 	Vertex s;
			// 	if (*itDirs != *itOrthDirs && *itDirs != *(++itOrthDirs)) continue;
			// 	if (computeMapInOneDirection(s, result, trackerNeighbors, itDirs, it->second).size() == 0) {
			// 		trace.info() << "returnin false" << endl;
			// 		return false;
			// 	}
			// }
			return true;
		}
	}
	return false;
}

template <typename MyDigitalSurface, typename Vertex>
vector<Vertex> extractPseudoBisector(MyDigitalSurface& boundary, const Vertex& start, const Vertex& end,  Viewer3D<>& viewer) {
	 typedef DiscreteExteriorCalculus<3, EigenLinearAlgebraBackend> DEC;
	 DEC calculus;
	 map<Vertex, vector<Vertex>> vertices;
	 //Creating the structure
	 for ( typename MyDigitalSurface::ConstIterator it = boundary.begin(), it_end = boundary.end();
		   it != it_end; ++it )
	 {
		 KSpace::SCells oneNeig = calculus.myKSpace.sLowerIncident(*it);
		 calculus.insertSCell(*it);
		 for(KSpace::SCells::ConstIterator itt = oneNeig.begin(), ittend = oneNeig.end(); itt != ittend; ++itt)
		 {
			 calculus.insertSCell(*itt);
			 KSpace::SCells oneNeig2 = calculus.myKSpace.sLowerIncident(*itt);
			 for(KSpace::SCells::ConstIterator ittt = oneNeig2.begin(), itttend = oneNeig2.end(); ittt != itttend; ++ittt) {
				 vertices[*it].push_back(*ittt);
				 calculus.insertSCell(*ittt);
			 }
		 }
	 }
  
	 //Setting a dirac on a 0-cell ==> 0-form
	 DEC::DualForm1 dirac(calculus);
	 DEC::Index ind = calculus.getSCellIndex(start);
	 dirac.myContainer( ind ) = 1;
  
	 //Laplace operator
	 DEC::DualDerivative1 dp1 = calculus.derivative<1,DUAL>();
	 DEC::DualHodge2 hodg2 = calculus.dualHodge<2>();
	 DEC::PrimalDerivative1 d1 = calculus.derivative<1, PRIMAL>();
	 DEC::PrimalHodge2 phodg2 = calculus.primalHodge<2>();

	 //Diffusion-like operator
	 DEC::DualIdentity1 laplace=   calculus.identity<1, DUAL>() - phodg2*d1*hodg2*dp1 ; //calculus.primalLaplace() ;
	 DEC::DualIdentity1 lap2 = laplace*laplace;
	 DEC::DualIdentity1 lap4 = lap2*lap2;
	 DEC::DualIdentity1 lap8 = lap4*lap4;
	 DEC::DualIdentity1 lap16 = lap8*lap8;
  
	 //Solver
	 typedef EigenLinearAlgebraBackend::SolverSimplicialLLT LinearAlgebra;
	 typedef DiscreteExteriorCalculusSolver<DEC, LinearAlgebra, 1, DUAL, 1, DUAL> Solver;
	 Solver solver;
	 solver.compute(lap2);
	 DEC::DualForm1 result = solver.solve(dirac);
	 trace.info() << result << endl;
	 
	 DEC::PrimalVectorField vf = calculus.sharp(hodg2 * dp1 * result);
	 DEC::PrimalVectorField normvf = vf.normalized();
	
//	 Display3DFactory<>::draw(viewer, normvf);
	 vector<Vertex> path;
	 
   
	 /*for (typename DEC::Index index=0; index<result.myCalculus->kFormLength(1, DUAL); index++)
	   {
		 Vertex vertex = result.getSCell(index);		 
		 const bool& flipped = result.myCalculus->isSCellFlipped(vertex);
		 if (flipped) vertex = result.myCalculus->myKSpace.sOpp(vertex);
		 vector<Vertex> correspondingVertices = vertices.at(vertex);
		 if (vertex == start) continue;
		 
		 vector<RealPoint> values;
		 for (int i =0; i < 4; i++) {
			 RealPoint arrow = normvf.getArrow(normvf.myCalculus->getSCellIndex(correspondingVertices[i]));
			 values.push_back(arrow);
			 if (vertex == end) trace.info() << arrow << endl;
		 }
		 for (int i = 0; i < values.size() - 1; i++) {
			 for (int j = i+1; j < values.size(); j++) {
				 if (allSigns(values[i], values[j])) {
					 path.push_back(vertex);
				 }
			 }
		 }
		 }*/
	
	 for (typename DEC::Index index=0; index<result.myCalculus->kFormLength(1, DUAL); index++)
	 {
		 int cpt = 0;
		 Vertex vertex = result.getSCell(index);
		 const bool& flipped = result.myCalculus->isSCellFlipped(vertex);
		 if (flipped) vertex = result.myCalculus->myKSpace.sOpp(vertex);
		 DEC::DualForm1::Scalar currentScalar = result.myContainer(index);
		 if (flipped) currentScalar = -currentScalar;
		 vector<Vertex> neighbors;
		 back_insert_iterator<vector<Vertex>> iter(neighbors);
		 boundary.writeNeighbors(iter, vertex);
		 for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			 DEC::Index index = result.myCalculus->getSCellIndex(*it);
			 const bool& flipped = result.myCalculus->isSCellFlipped(*it);
			 DEC::DualForm1::Scalar otherScalar = result.myContainer(index);
			 if (flipped) otherScalar = -otherScalar;
			 if (currentScalar < otherScalar) cpt++;
		  }
		  if (cpt == 4) {
			  viewer << CustomColors3D(Color::Cyan, Color::Cyan) << vertex;
			 path.push_back(vertex);
		  }
		 
	 }

	 draw(viewer, result, 0, 0);
	 /*for (auto it = path.begin(), ite = path.end(); it != ite; ++it) {
		 int cpt = 0;
		 vector<Vertex> neighbors;
		 back_insert_iterator<vector<Vertex>> iter(neighbors);
		 boundary.writeNeighbors(iter, *it);
		 for (auto itn = neighbors.begin(), itne = neighbors.end(); itn != itne; ++itn) {
			 if (find(path.begin(), path.end(), *itn) != path.end()) {
				 cpt++;
			 }
		 }
		 if (cpt == 2) { //needs to be a curve, two neighbors
			 viewer << CustomColors3D(Color::Cyan, Color::Cyan) << *it;
			 path2.push_back(*it);
		 }
		 }*/

		
	 
		   
	 
	 return path;
}

template <typename Visitor, typename Vertex>
vector<Vertex> extractGeodesicLoop(Visitor& visitor, const vector<Vertex>& pseudoBisector, const Vertex& bel) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	map<Vertex, Vertex> aMapPrevious;
	vector<Vertex> thePath;
	
	while (!visitor.finished()) {
		node = visitor.current();
		vector<Vertex> neighbors;
		back_insert_iterator<vector<Vertex>> iter(neighbors);
		visitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			auto itMapExisting = aMapPrevious.find(*it);
			if (itMapExisting == aMapPrevious.end()) {
				aMapPrevious[*it] = node.first;
			}
			else {
				if (find(pseudoBisector.begin(), pseudoBisector.end(), *it) == pseudoBisector.end())
					continue;
				Vertex tmp = node.first;
			    Vertex tmp2 = aMapPrevious[*it];
				set<Vertex> path1, path2;
				vector<Vertex> vpath1, vpath2;
				Vertex nearBel, otherNearBel;
				while (tmp != bel) {
					path1.insert(tmp);
					vpath1.push_back(tmp);
					tmp = aMapPrevious[tmp];
					if (aMapPrevious[tmp] == bel) nearBel = tmp;
				}
				while (tmp2 != bel) {
					path2.insert(tmp2);
					vpath2.push_back(tmp2);
					tmp2 = aMapPrevious[tmp2];
					if (aMapPrevious[tmp2] == bel) otherNearBel = tmp2;
				}
				set<Vertex> intersect;
				set_intersection(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(intersect, intersect.begin()));
				if (intersect.size() != 0) continue;
				vector<Vertex> correspondingPath;
				
				
				std::copy (vpath1.rbegin(), vpath1.rend(), std::back_inserter(correspondingPath));
				correspondingPath.push_back(bel);
				std::copy (vpath2.begin(), vpath2.end(), std::back_inserter(correspondingPath));
				correspondingPath.push_back(*it);
				
			    thePath = correspondingPath;
				visitor.terminate();
			}
		}
		if (!visitor.finished())
			visitor.expand();
	}
	return thePath;
}

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}


/**
 * Returns a surfel associated with a Z3 point, using a container mapping surfels to points
 */
template <typename Surfel, typename Point>
Surfel convertToSurfel(const Point& point, const map<Surfel, Point> & surfelToPoint) {
	Surfel surfel;
	for (auto it = surfelToPoint.begin(), ite = surfelToPoint.end(); it != ite; ++it) {
		if (it->second == point) {
			return it->first;
		}
	}
	return surfel;
}

int main( int argc, char** argv )
{
	typedef PointVector<3, double> Vector;
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef ImageSelector<Z3i::Domain, float>::Type FloatImage;
	typedef functors::IntervalForegroundPredicate<Image> Binarizer;
	typedef functors::EqualPointPredicate<Point> EqualPredicate;
	typedef functors::NotPointPredicate<Binarizer> NotPredicate;
	typedef Z3i::KSpace KSpace;
	typedef LightImplicitDigitalSurface<KSpace, Z3i::DigitalSet >   MyDigitalSurfaceContainer;
	typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
	typedef SurfacePoint<Z3i::Point, SCell> SurfacePoint;
	
	typedef BreadthFirstVisitor<MyDigitalSurface> MyBreadthFirstVisitor;
	typedef MyBreadthFirstVisitor::Node MyNode;
	typedef MyBreadthFirstVisitor::Size MySize;
	typedef Path<SCell> Path;
	
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file")
		("skeleton,s", po::value<std::string>(), "vol file (skeleton)")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
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
	string inputSkeletonName = vm["skeleton"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();



	Image img  = VolReader<Image>::importVol( inputFilename );
	Image skeleton = VolReader<Image>::importVol( inputSkeletonName );
	Z3i::Domain domain = img.domain();
	KSpace ks;
	ks.init( domain.lowerBound(), 
			 domain.upperBound(), true );

	QApplication app(argc, argv);
	Viewer3D<> viewer(ks);

	Cell dummy;
	viewer << SetMode3D( dummy.className(), "Basic" );
	viewer.show();
	
	Z3i::DigitalSet set3d (domain);
	SetFromImage<Z3i::DigitalSet>::append<Image> (set3d, img, 
												  thresholdMin, thresholdMax);
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 10000000 );
	typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
	MySurfelAdjacency surfAdj( true );
	MyDigitalSurfaceContainer* ptrSurfContainer = 
		new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
	MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
	set<SCell> surfelSet;
	map<SCell, Point> surfaceVoxelSet;
	set<Point> surfacePointSet;
	for (auto it = digSurf.begin(), itE = digSurf.end(); it!=itE; ++it) {
		surfelSet.insert(*it);
	}
	
	SurfaceUtils::computeSurfelWeight<SurfacePoint>(ks, set3d, surfelSet, surfacePointSet, surfaceVoxelSet);


	Z3i::DigitalSet setskel (skeleton.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image> (setskel, skeleton, 
												  thresholdMin, thresholdMax);
	Color color(100,100,140,128);
    int number = 1;
	int i = 0;
	trace.beginBlock("Computing distance map");
	while (i < number) {
		Domain domainskel = skeleton.domain();
		SCell surfel, endSurfel;
		do {
			auto it = select_randomly(setskel.begin(), setskel.end());
			Point center = *it;
			Point point = pointCorrespondingToMinDistance<Point>(center, surfacePointSet);
			Point projected = computeProjection(point, center, set3d);
	
		    surfel = convertToSurfel(point, surfaceVoxelSet);
			endSurfel = convertToSurfel(projected, surfaceVoxelSet);
		} while (surfel == SCell() || endSurfel == SCell());
		
		
		trace.info() << surfel << " " << endSurfel << endl;
		
		//	surfel = {{-36, 1, 5}, true};
		//endSurfel = {{42, 5, 5}, false};
		viewer << CustomColors3D(Color::Green, Color::Green) << surfel;
		viewer << CustomColors3D(Color::Yellow, Color::Yellow) << endSurfel;
		vector<SCell> pseudoBisector = extractPseudoBisector(digSurf, surfel, endSurfel, viewer);

		for (auto it = pseudoBisector.begin(), ite = pseudoBisector.end(); it != ite; ++it) {
			viewer << CustomColors3D(Color::Cyan, Color::Cyan) << *it;
		}
		MyBreadthFirstVisitor visitor( digSurf, surfel );
		vector<SCell> path = extractGeodesicLoop(visitor, pseudoBisector, surfel);

		for (auto it = path.begin(), ite = path.end(); it != ite; ++it) {
			viewer << CustomColors3D(Color::Red, Color::Red) << *it;
		}
	    i++;
	}
	trace.endBlock();

	trace.beginBlock("Displaying");
	for (auto it = set3d.begin(), ite = set3d.end(); it != ite; ++it) {
		//		viewer << CustomColors3D(color, color) << *it;
	}
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
	trace.endBlock();
	
	return 0;
}
