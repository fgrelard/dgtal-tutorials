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
#include "DGtal/dec/DiscreteExteriorCalculusFactory.h"

#include "surface/SurfaceUtils.h"
#include "surface/SurfacePoint.h"
#include "geometry/Distance.h"
#include "shapes/Ball.h"
#include "geometry/PointUtil.h"
#include "geometry/SliceUtils.h"
#include "shapes/Path.h"

using namespace std;
using namespace DGtal;
using namespace Z3i;
namespace po = boost::program_options;


template <typename Vector>
bool allSigns(const Vector& v1, const Vector& v2, unsigned int dir) {
	double max = -1, max2 = -1;
	int index = -1, index2 = -2;
	for (int i = 0; i < 3; i++) {
		if (abs(v1[i]) > max) {
			index = i;
			max = abs(v1[i]);
		}
	}
	for (int i = 0; i < 3; i++) {
		if (abs(v2[i]) > max2) {
			index2 = i;
			max2 = abs(v2[i]);
		}
	}
    return (index == index2 &&  v1[index] * v2[index] < 0);
/*	double x = v1[0] * v2[0];
	double y = v1[1] * v2[1];
	double z = v1[2] * v2[2];
	return (x <= 0 && y <= 0 && z <= 0);*/
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
			if (PointUtil::areNeighbors(*it, *secondIt) && *secondIt != *it) {
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
                            
    for (typename TCalculus::Index index=0; index<kform.length(); index++)
    {
        SCell displayed_cell =  kform.getSCell(index);
        Scalar displayed_value = kform.myContainer(index);
		
        if (!std::isfinite(log10(displayed_value))) displayed_value = -displayed_value;

        if (std::isfinite(displayed_value)) {
			display << DGtal::CustomColors3D(DGtal::Color::Black, colormap(abs(log10(displayed_value))/denominator));
		}
        else
            display << DGtal::CustomColors3D(DGtal::Color::Black, DGtal::Color::White);

        display << displayed_cell;
    }
}

template <typename Vertex, typename DEC, DGtal::Order order, DGtal::Duality duality>
map<Vertex, double> constructLevelSetMap(const DGtal::KForm<DEC, order, duality>& kform, double minValue, double maxValue) {
	typedef typename DEC::Scalar Scalar;
    typedef typename DEC::SCell SCell;
    typedef typename DEC::Index Index;

	map<Vertex, double> mapLevelSet;
	
	double denominator = 0.0;
	bool first = true;
	for (Index index=0; index<kform.myContainer.rows(); index++)
	{
		const Scalar value = kform.myContainer(index);
		if (!std::isfinite(abs(log10(value)))) continue;
		if (first || denominator < abs(log10(value))) denominator = abs(log10(value));
		first = false;
	}

	for (Index index=0; index<kform.myCalculus->kFormLength(order, duality); index++)
    {
        const SCell& cell = kform.myCalculus->getSCell(order, duality, index);
        ASSERT(kform.myCalculus->myKSpace.sSign(cell) == DEC::KSpace::POS);

        const bool& flipped = kform.myCalculus->isSCellFlipped(cell);
        SCell displayed_cell = cell;
        if (flipped) displayed_cell = kform.myCalculus->myKSpace.sOpp(cell);
		
		Scalar displayed_value = kform.myContainer(index);
        if (!std::isfinite(log10(displayed_value))) displayed_value = -displayed_value;

		//level set
		displayed_value = roundf(abs(log10(displayed_value))/denominator * 10)/10;
		
		mapLevelSet[displayed_cell] = displayed_value;
		
	}
	return mapLevelSet;
	
}

template <typename Vertex, typename Result, typename Tracker, typename DirIterator>
map<Vertex, double> computeMapInOneDirection(Vertex& s,  const Result& result, Tracker tracker, DirIterator& itDirs,double currentScalar ) {
	map<Vertex, double> aMapVertexToScalar;
	if ( tracker->adjacent( s, *itDirs, true )) {
		typename Result::Calculus::Index index = result.myCalculus->getSCellIndex(s);
		const bool& flipped = result.myCalculus->isSCellFlipped(s);
		typename Result::Calculus::DualForm1::Scalar otherScalar = result.myContainer(index);
		if (flipped) otherScalar = -otherScalar;
		if (currentScalar < otherScalar) aMapVertexToScalar[s] = otherScalar;
	}
	if ( tracker->adjacent( s, *itDirs, false )) {
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
	auto itDirs  = result.myCalculus->myKSpace.sDirs(vertex);
	auto tracker = boundary.container().newTracker( vertex );
	tracker->move(vertex);
	Vertex s;
	for (; itDirs != 0; ++itDirs) {
	    map<Vertex, double> aMapVertexToScalar = computeMapInOneDirection(s, result, tracker, itDirs, currentScalar);
		if (aMapVertexToScalar.size() == 2) {
		
			Vertex first = (aMapVertexToScalar.begin())->first;
			auto itDirs =  result.myCalculus->myKSpace.sDirs(first);
			auto tracker = boundary.container().newTracker(first);
			tracker->move(first);
			Vertex s;
			bool passing = false;
			for (; itDirs != 0; ++itDirs) {
				map<Vertex, double> amap = computeMapInOneDirection(s, result, tracker, itDirs, currentScalar);
				if (amap.size() == 2) {
					passing = true;
				}
			}
			if (!passing) continue;
			Vertex second = (++aMapVertexToScalar.begin())->first;
			itDirs =  result.myCalculus->myKSpace.sDirs(second);
			tracker = boundary.container().newTracker(second);
			tracker->move(second);
			s = Vertex();
			passing = false;
			for (; itDirs != 0; ++itDirs) {
				map<Vertex, double> amap = computeMapInOneDirection(s, result, tracker, itDirs, currentScalar);
				if (amap.size() == 2) {
					passing = true;
				}
			}
			if (passing) {
				return true;
				Vertex first = (aMapVertexToScalar.begin())->first;
				Vertex second = (++aMapVertexToScalar.begin())->first;
				vector<Vertex> neighFirst;
				back_insert_iterator<vector<Vertex>> itfirst(neighFirst);
				boundary.writeNeighbors(itfirst, first);

			
				vector<Vertex> neighSecond;
				back_insert_iterator<vector<Vertex>> itsecond(neighSecond);
				boundary.writeNeighbors(itsecond, second);

				for (auto it = neighFirst.begin(), ite = neighFirst.end(); it != ite; ++it) {
					for (auto its = neighSecond.begin(), itse = neighSecond.end(); its != itse; ++its) {
						if (*it == *its && *it != vertex)
							return false;
					}
				}
				return true;
			}
		}
	}
	return false;
}

template <typename MyDigitalSurface, typename Vertex>
vector<Vertex> extractPseudoBisector(MyDigitalSurface& boundary, const Vertex& start, const Vertex& end,  Viewer3D<>& viewer) {
	typedef DGtal::DiscreteExteriorCalculusFactory<DGtal::EigenLinearAlgebraBackend> CalculusFactory;
	typedef DiscreteExteriorCalculus<2, 3, EigenLinearAlgebraBackend> DEC;
	const DEC calculus = CalculusFactory::createFromNSCells<2>(boundary.begin(), boundary.end());
	// for ( typename MyDigitalSurface::ConstIterator it = boundary.begin(), it_end = boundary.end();
	// 	   it != it_end; ++it )
	//  {
	// 	 KSpace::SCells oneNeig = calculus.myKSpace.sLowerIncident(*it);
	// 	 calculus.insertSCell(*it);
	// 	 for(KSpace::SCells::ConstIterator itt = oneNeig.begin(), ittend = oneNeig.end(); itt != ittend; ++itt)
	// 	 {
	// 		 calculus.insertSCell(*itt);
	// 		 KSpace::SCells oneNeig2 = calculus.myKSpace.sLowerIncident(*itt);
	// 		 for(KSpace::SCells::ConstIterator ittt = oneNeig2.begin(), itttend = oneNeig2.end(); ittt != itttend; ++ittt) {
	// 			 calculus.insertSCell(*ittt);
	// 		 }
	// 	 }
	//  }
	 trace.info() << calculus << endl;
	 //Setting a dirac on a 0-cell ==> 0-form
	 DEC::DualForm0 dirac(calculus);
	 DEC::Index ind = calculus.getCellIndex(calculus.myKSpace.unsigns(start));
	 dirac.myContainer( ind ) = 1;
  
	 //Laplace operator
	 const DEC::DualIdentity0 laplace = calculus.laplace<DUAL>();
	 const DEC::DualDerivative0 d0p = calculus.derivative<0, DUAL>();
	 const DEC::DualDerivative1 d1p = calculus.derivative<1, DUAL>();
	 const DEC::DualHodge1 hodg1p = calculus.hodge<1, DUAL>();
	 typedef DGtal::EigenLinearAlgebraBackend::SolverSparseLU LinearAlgebraSolver;
	 typedef DGtal::DiscreteExteriorCalculusSolver<DEC, LinearAlgebraSolver, 0, DUAL, 0, DUAL> PoissonSolver;

	 PoissonSolver solver;
	 solver.compute(laplace);
	 DEC::DualForm0 result = solver.solve(dirac);
	 trace.info() << result << endl;
	 DEC::PrimalForm1 oneFormResult = hodg1p * d0p * result;
	 DEC::PrimalVectorField vf = calculus.sharp(oneFormResult);
	 DEC::PrimalVectorField normvf = vf.normalized();
//	 Display3DFactory<Space, KSpace>::draw(viewer, oneFormResult);
	 Display3DFactory<Space, KSpace>::draw(viewer, normvf);
	 vector<Vertex> path;
	 
   
	 // for (typename DEC::Index index=0; index<result.myCalculus->kFormLength(1, DUAL); index++)
	 // {
	 // 	 Vertex vertex = result.getSCell(index);		 
	 // 	 const bool& flipped = result.myCalculus->isSCellFlipped(vertex);
	 // 	 if (flipped) vertex = result.myCalculus->myKSpace.sOpp(vertex);
	 // 	 vector<Vertex> correspondingVertices = vertices.at(vertex);
	 // 	 if (vertex == start) continue;

	 // 	 unsigned int dir = result.myCalculus->myKSpace.sOrthDir(vertex);
		 
	 // 	 vector<RealPoint> values;
	 // 	 for (int i =0; i < 4; i++) {
	 // 		 RealPoint arrow = normvf.getArrow(normvf.myCalculus->getSCellIndex(correspondingVertices[i]));
	 // 		 values.push_back(arrow);
	 // 		 if (vertex == end) trace.info() << arrow << endl;
	 // 	 }
	 // 	 if (allSigns(values[0], values[3], dir) && allSigns(values[1], values[2], dir)) {
	 // 		 path.push_back(vertex);
	 // 	 }			 
	 // }
	 
	 int size = result.myCalculus->kFormLength(0, DUAL);
	 for (typename DEC::Index index=0; index<size; index++)
	 {
		 int cpt = 0;
	 	 Vertex vertex = result.getSCell(index);
	 	 DEC::DualForm0::Scalar currentScalar = result.myContainer(index);
		 vector<Vertex> neighbors;
	 	 back_insert_iterator<vector<Vertex>> iter(neighbors);
	 	 boundary.writeNeighbors(iter, vertex);
		 for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			 DEC::Index index = result.myCalculus->getCellIndex(calculus.myKSpace.unsigns(*it));
	 		 DEC::DualForm0::Scalar otherScalar = result.myContainer(index);
	 		 if (currentScalar < otherScalar) cpt++;
	 	 }
	 	 if (cpt == 4)
			 path.push_back(vertex);
		 /* if ( belongsToBisector(vertex, result, boundary, currentScalar)) {
	 		 path.push_back(vertex);
	 		 }*/
		 
	 }
	 // /*vector<Vertex> path2;
	 // for (auto it = path.begin(); it != path.end(); ++it) {
	 // 	 int cpt = 0;
	 // 	 vector<Vertex> neighbors;
	 // 	 back_insert_iterator<vector<Vertex>> iter(neighbors);
	 // 	 boundary.writeNeighbors(iter, *it);
	 // 	 for (auto itn = neighbors.begin(), itne = neighbors.end(); itn != itne; ++itn) {
	 // 		 if (find(path.begin(), path.end(), *itn) != path.end()) {
	 // 			 cpt++;
	 // 		 }
	 // 	 }
	 // 	 if (cpt == 1) { //needs to be a curve, at least one neighbor for it to be a pseudo bisector
	 // 		 path2.push_back(*it);
	 // 	 }
	 // }*/
		 
	 // draw(viewer, result, 0, 0);
		
 
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

	DigitalSet surfacePointDigitalSet(domain);
	surfacePointDigitalSet.insert(surfacePointSet.begin(), surfacePointSet.end());
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
			Point projected = computeProjection(point, center, surfacePointDigitalSet);
	
		    surfel = convertToSurfel(point, surfaceVoxelSet);
			endSurfel = convertToSurfel(projected, surfaceVoxelSet);
		} while (surfel == SCell() || endSurfel == SCell());
		
		trace.info() << surfel << " " << endSurfel << endl;
		
		//surfel = {{-36, 1, 5}, true};
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
		//	viewer << CustomColors3D(color, color) << *it;
	}
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
	trace.endBlock();
	
	return 0;
}
