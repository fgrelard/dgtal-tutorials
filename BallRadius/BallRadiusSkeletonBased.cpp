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
	return (x < 0 && y < 0 && z < 0);
}

template <typename Path, typename Vertex>
bool
computeFlowOnPath(const Path& path, const Point & center, const Vertex& bel, Viewer3D<>& viewer) {
	typedef DiscreteExteriorCalculus<3, EigenLinearAlgebraBackend> Calculus;
	Calculus calculus;
	
	double mean = 0;
	const Calculus::DualHodge1 h2 = calculus.dualHodge<1>();
	vector<Vertex> points = path.myPath;
		
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		auto edge = calculus.myKSpace.sLowerIncident(*it);
	   	calculus.insertSCell(*it);
	    for (auto eit = edge.begin(), ee = edge.end(); eit != ee; ++eit) {
			auto point = calculus.myKSpace.sLowerIncident(*eit);
		   	calculus.insertSCell(*eit);
			for (auto pit = point.begin(), pite = point.end(); pit != pite; ++pit) 
				calculus.insertSCell(*pit);
		}
	}

	vector<Vertex> vertices;
	auto edge = calculus.myKSpace.sLowerIncident(bel);
	for (auto eit = edge.begin(), ee = edge.end(); eit != ee; ++eit) {
		auto point = calculus.myKSpace.sLowerIncident(*eit);
		for (auto pit = point.begin(), pe = point.end(); pit != pe; ++pit) {
			vertices.push_back(*pit);
		}
	}

	Calculus::DualForm1 dirac(calculus);
	Calculus::DualDerivative1 dp1 = calculus.derivative<1,DUAL>();
	Calculus::DualHodge2 hodg2 = calculus.dualHodge<2>();
	Calculus::PrimalDerivative1 d1 = calculus.derivative<1, PRIMAL>();
	Calculus::PrimalHodge2 phodg2 = calculus.primalHodge<2>();
 
	//For gradient
	const Calculus::PrimalDerivative0 d0 = calculus.derivative<0, PRIMAL>();
	//Diffusion-like operator
    Calculus::DualIdentity1 laplace = calculus.identity<1, DUAL>() - phodg2*d1*hodg2*dp1 ; //calculus.primalLaplace() ;
	
	PointVector<3, double> p;
	if (dirac.myContainer.size() > 1) {
		Calculus::Index ind =  dirac.myCalculus->getSCellIndex(bel);
		dirac.myContainer(ind) = 1;

		typedef EigenLinearAlgebraBackend::SolverSimplicialLLT LinearAlgebraSolver;
		typedef DiscreteExteriorCalculusSolver<Calculus, LinearAlgebraSolver, 1, DUAL, 1, DUAL> Solver;
		Solver solver;
		solver.compute(laplace);
		Calculus::DualForm1 solved_solution = solver.solve(dirac);
	    int size = solved_solution.myContainer.innerSize();
		for (int i  = 0; i < size; i++) {
			mean += solved_solution.myContainer[i];
		}
		mean /= size;
	    VectorField<Calculus, PRIMAL> vf = calculus.sharp(hodg2 * dp1 * solved_solution);
		VectorField<Calculus, PRIMAL> normvf = vf.normalized();
//		Display3DFactory<Space, KSpace>::draw(viewer, normvf);

		vector<PointVector<3, double>> cross;
		for (auto it = vertices.begin(), ite = vertices.end(); it != ite; ++it) {
			Calculus::Index ind = normvf.myCalculus->getSCellIndex(*it);
			DGtal::Z3i::RealPoint origin = normvf.myCalculus->myKSpace.sKCoords(*it)/2.;
			origin -= {.5, .5, .5};
			PointVector<3, double> toNext = normvf.getArrow(ind);
			PointVector<3, double> toCenter = PointVector<3, double>((center - origin)).getNormalized();
			PointVector<3, double> crossProduct = toNext.crossProduct(toCenter);
			cross.push_back(crossProduct);
		}
		for (auto it = cross.begin(), ite = cross.end(); it != ite; ++it) {
			for (auto otherIt = cross.begin(), otherite = cross.end(); otherIt != otherite; ++otherIt) {
				if (allSigns(*it, *otherIt))
					return true;
			}
		}
		/*for (typename Calculus::Index index=0; index<normvf.myCalculus->kFormLength(0, PRIMAL); index++)
		{
			p+= normvf.getArrow(index);
			const typename Calculus::SCell& cell = normvf.myCalculus->getSCell(0, PRIMAL, index);
			DGtal::Z3i::RealPoint origin = normvf.myCalculus->myKSpace.sKCoords(cell)/2.;
			origin -= {.5, .5, .5};
			PointVector<3, double> toNext = normvf.getArrow(index);
			PointVector<3, double> toCenter = PointVector<3, double>((center - origin)).getNormalized();
			PointVector<3, double> crossProduct = toNext.crossProduct(toCenter);
			if (crossProduct != PointVector<3, double>::zero) {
				viewer.addCone(origin, origin+crossProduct);
			}
			}*/
	}
	return false;
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



template <typename Path, typename Visitor, typename Vertex>
vector<Path> computeSystemOfLoops(Visitor& visitor, const Vertex& bel) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	map<Vertex, Vertex> aMapPrevious;
	vector<Path> systemOfLoops;
	
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
				
				Path path(correspondingPath, bel, *it, make_pair(nearBel, otherNearBel));
				systemOfLoops.push_back(path);
			}
		}
		if (!visitor.finished())
			visitor.expand();
	}
	return systemOfLoops;
}

template <typename Path>
Path computeMinimumGradientPath(const vector<Path>& paths, const Point& center, Viewer3D<>& viewer) {
	double minimum = std::numeric_limits<double>::max();
	Path path;
	for (auto it = paths.begin(), ite = paths.end(); it != ite; ++it) {	
		bool isOk = computeFlowOnPath(*it, center, it->myBel, viewer);
		if (isOk)
			return *it;
	}
	return path;
}

template <typename Path>
vector<Path> selectGeodesicLoops(const vector<Path>& systemOfLoops, const Path& path) {
	vector<Path> geodesicLoops;
	for (auto it = systemOfLoops.begin(), ite = systemOfLoops.end();
		 it != ite; ++it) {
		for (auto itPath = path.myPath.begin(), itPathE = path.myPath.end(); itPath != itPathE; ++itPath) {
			if (find(it->myPath.begin(), it->myPath.end(), *itPath) == it->myPath.end())
				continue;
		}
		geodesicLoops.push_back(*it);
	}
	return geodesicLoops;
	
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
bool sameSign(const vector<Point>& points) {
	int negcptX = 0, negcptY = 0, negcptZ = 0;
	int poscptX = 0, poscptY = 0, poscptZ = 0;
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		negcptX = ((*it)[0] < 0 ) ? negcptX + 1 : negcptX;
		negcptY = ((*it)[1] < 0 ) ? negcptY + 1 : negcptY;
		negcptZ = ((*it)[2] < 0 ) ? negcptZ + 1 : negcptZ;

		poscptX = ((*it)[0] > 0 ) ? poscptX + 1 : poscptX;
		poscptY = ((*it)[1] > 0 ) ? poscptY + 1 : poscptY;
		poscptZ = ((*it)[2] > 0 ) ? poscptZ + 1 : poscptZ;
	}
	bool x = (negcptX == points.size() || poscptX == points.size());
	bool y = (negcptY == points.size() || poscptY == points.size());
	bool z = (negcptZ == points.size() || poscptZ == points.size());
	return ( (x && y) || (x && z) || (y && z));
}



template <typename Point>
bool consistentCrossProductsAlong(const vector<Point>& path, const Point& center) {
	map<Point, Point> pointToDirection = computeDSSOnPath(path);
	vector<Point> crossProducts;
	for (auto it = pointToDirection.begin(), ite = pointToDirection.end(); it != ite; ++it) {
		Point vectorToCenter = center - it->first;
		Point vectorToNext = it->second;
		Point crossProduct = vectorToNext.crossProduct(vectorToCenter);
		crossProducts.push_back(crossProduct);
	}
	return sameSign(crossProducts);
}

template <typename Surfel, typename Point>
vector<Point> selectGeodesicLoops(const vector<Path<Surfel>>& systemOfLoops, const map<Surfel, Point>& surfelToPoint, const Point& center) {
	vector<Point> geodesicLoops;
	for (auto it = systemOfLoops.begin(), ite = systemOfLoops.end(); it != ite; ++it) {
		vector<Point> pathPoints;
		for (auto itPath = it->begin(), itPathE = it->end(); itPath != itPathE; ++itPath) {
			pathPoints.push_back(surfelToPoint.at(*itPath));
		}
		bool consistentCrossP = consistentCrossProductsAlong(pathPoints, center);
		if (consistentCrossP)
			std::copy (pathPoints.begin(), pathPoints.end(), std::back_inserter(geodesicLoops));
	}
	return geodesicLoops;
}

template <typename Domain>
bool isSameProjectedDomain(const Domain& domain, const typename Domain::Point& center) {
	typedef typename Domain::Point Point;
	Point mini = domain.upperBound(), maxi = domain.lowerBound();
	for (typename Domain::Iterator it = domain.begin(), ite = domain.end(); it!= ite; ++it) {
		Point toCenter = center - *it;
		Point fromCenter = center + toCenter;
		if (fromCenter < mini) mini = fromCenter;
		if (fromCenter > maxi) maxi = fromCenter;
	}
	return mini == domain.lowerBound() && maxi == domain.upperBound();
}

template <typename Domain, typename Vector, typename Point, typename Visitor, typename Vertex>
vector<Point> createShortestPath(Visitor& visitor, const Vertex& bel, const Point& center, const map<Vertex, Point>& surfelToPoint) {
	typedef typename Visitor::Node MyNode;
	
	MyNode node;
	map<Vertex, Vertex> aMapPrevious;
	vector<Point> thePath;
	
	while (!visitor.finished()) {
		if (node.second == 100) break;
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
			    Vertex tmp = node.first;
			    Vertex tmp2 = aMapPrevious[*it];
				vector<Point> correspondingPoints;
				set<Vertex> path1, path2;

				while (tmp != bel) {
					path1.insert(tmp);
					tmp = aMapPrevious[tmp];
				}
				while (tmp2 != bel) {
					path2.insert(tmp2);
					tmp2 = aMapPrevious[tmp2];
				}

				//creating corresponding path with voxels
				correspondingPoints.push_back(surfelToPoint.at(bel));
				for (auto it = path1.rbegin(), ite = path1.rend(); it != ite; ++it) {
					correspondingPoints.push_back(surfelToPoint.at(*it));
				}
				correspondingPoints.push_back(surfelToPoint.at(*it));
				for (auto it = path2.begin(), ite = path2.end(); it != ite; ++it) {
					correspondingPoints.push_back(surfelToPoint.at(*it));
				}

				//checking if the intersect is null
				set<Vertex> intersect;
				set_intersection(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(intersect, intersect.begin()));
				if (intersect.size() != 0) continue;
			    Vector normal = SliceUtils::computeNormalFromCovarianceMatrix<Vector>(correspondingPoints);
				Domain domain = PointUtil::computeBoundingBox<Domain>(correspondingPoints);
				bool isProjected = isSameProjectedDomain(domain, center);
				
				//checks if the center point belongs to the plane
				//double eq = normal[0] * center[0] +
				//	normal[1] * center[1] +
				//	normal[2] * center[2];
				if (isProjected && normal != Vector::zero) {
					visitor.terminate();		  
					thePath = correspondingPoints;
				}
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
    int number = 10;
	int i = 0;
	trace.beginBlock("Computing distance map");
	while (i < number) {
		Domain domainskel = skeleton.domain();
		auto it = select_randomly(setskel.begin(), setskel.end());	 
		Point point = pointCorrespondingToMinDistance<Point>(*it, surfacePointSet);
		SCell surfel = convertToSurfel(point, surfaceVoxelSet);
		trace.info() << surfel << endl;
		viewer << CustomColors3D(Color::Green, Color::Green) << surfel;
		if (surfel != SCell()) {
			MyBreadthFirstVisitor visitor( digSurf, surfel );
			/*vector<Path> systemOfLoops = computeSystemOfLoops<Path>(visitor, surfel);
			  Path thePath = computeMinimumGradientPath(systemOfLoops, *it, viewer);
			  vector<Path> geodesicLoops = selectGeodesicLoops(systemOfLoops, thePath);
			  vector<SCell> scellsInPath = thePath.myPath;
			  for (auto it = scellsInPath.begin(), ite = scellsInPath.end(); it != ite; ++it) {
			  viewer << CustomColors3D(Color::Red, Color::Red) << *it;
			  }
			  for (auto it = geodesicLoops.begin(), ite = geodesicLoops.end(); it != ite; ++it) {
			  trace.info() << it->myPath.size() << endl;
			  for (auto itpath = it->begin(), itpathe = it->end(); itpath != itpathe; ++itpath) {
			  viewer << CustomColors3D(Color::Red, Color::Red) << *itpath;
			  }
			  }*/
			vector<Point> path = createShortestPath<Domain, Vector>(visitor, surfel, *it, surfaceVoxelSet);
			if (path.size() == 0) continue;
			for (auto it = path.begin(), ite = path.end(); it != ite; ++it) {
				viewer << CustomColors3D(Color::Red, Color::Red) << *it;
			}
			i++;
		}
	}
	trace.endBlock();

	trace.beginBlock("Displaying");
	for (auto it = set3d.begin(), ite = set3d.end(); it != ite; ++it) {
		viewer << CustomColors3D(color, color) << *it;
	}
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
	trace.endBlock();
	
	return 0;
}
