#include <iostream>
#include <thread>
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
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/TangentUtils.h"
#include "geometry/SliceUtils.h"
#include "geometry/Pencil.h"
#include "geometry/MSTTangent.h"
#include "geometry/PointUtil.h"
#include "geometry/WeightedPoint.h"
#include "surface/SurfaceUtils.h"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "Statistics.h"
#include "geometry/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "shapes/Ball.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

template <typename VCM, typename Domain>
Z3i::DigitalSet computeBranchingPartsWithVCMFeature(const VCM& vcm, const Domain& domain, double threshold) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	Z3i::DigitalSet aSet(domain);
	for (P2EConstIterator  it = vcm.mapPoint2ChiVCM().begin(),
			 itE = vcm.mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = lambda[ 0 ] / ( lambda[ 0 ] + lambda[ 1 ] + lambda[ 2 ] );
		if (ratio > threshold)
			aSet.insert(it->first);
	}
	return aSet; 
}

Z3i::Point extractNearestNeighborInSetFromPoint(const Z3i::DigitalSet& aSet, const Z3i::RealPoint& aPoint) {
	double distanceMin = numeric_limits<double>::max();
	Z3i::Point toReturn;
	for (auto it = aSet.begin(), ite = aSet.end(); it != ite; ++it) {
		double distanceToPoint = sqrt(pow(((*it)[0] - aPoint[0]), 2) + pow(((*it)[1] - aPoint[1]), 2) + pow(((*it)[2] - aPoint[2]), 2));
		if (distanceToPoint < distanceMin) {
			distanceMin = distanceToPoint;
			toReturn = *it;
		}
	}
	return toReturn;
}

template <typename Domain>
Z3i::DigitalSet linkTwoPointsWithTangents(const Z3i::Point& start,
										  const Z3i::Point& end,
										  const Z3i::RealPoint& tStart,
										  const Z3i::RealPoint& tEnd,
										  const Domain& domain) {
	Z3i::DigitalSet aSet(domain);
	Z3i::RealPoint vector = tStart;
	Z3i::RealPoint current = start;
	double initialDistance = Z3i::l2Metric(start, end)/2;
	while (Z3i::l2Metric(current, end) > sqrt(3) && Z3i::l2Metric(current, end) <= Z3i::l2Metric(start, end)) {		
	    current += vector;
		aSet.insert(current);
		double distanceToStart = std::abs(start.norm()- current.norm()) - initialDistance;
		if (distanceToStart < 0)
			distanceToStart = 1;
		double distanceToEnd = -std::abs(end.norm()- current.norm()) +  initialDistance;
		double weightStart = 1/distanceToStart;
		double weightEnd;
		if (distanceToEnd < 0)
			weightEnd = 0;
		else 
			weightEnd = 1/distanceToEnd;
		vector = ((weightStart * tStart + weightEnd * tEnd) / (weightStart + weightEnd)).getNormalized();
	}
	return aSet;	
}


Z3i::DigitalSet computeTraversedPoints(const Z3i::DigitalSet& setVolume,
									   const Z3i::Point& point,
									   const Z3i::RealPoint& normal ) {
	int scalar = 1;
	Z3i::DigitalSet traversed(setVolume.domain());
	
	Z3i::Point projection = point;
	while (setVolume.find(projection) != setVolume.end()) {
		projection = point + normal * scalar;
		traversed.insert(projection);
		scalar++;	
	}
	return traversed;
}

Z3i::DigitalSet computeSubVolume(const Z3i::DigitalSet& setVolume,
								 const Z3i::Point& point,
								 const Z3i::RealPoint& normal,								 
								 double radius) {
	Z3i::DigitalSet subVolume(setVolume.domain());
	int scalar = 1;
	Z3i::RealPoint theNormal = normal;
	Z3i::Point newPoint = point;
	while (setVolume.find(newPoint) != setVolume.end()) {
		Ball<Z3i::Point> ball(newPoint, radius);
		newPoint = point + theNormal * scalar;
		std::vector<Z3i::Point> pointsInBall = ball.pointsInBall();
		for (auto it = pointsInBall.begin(), ite = pointsInBall.end(); it !=ite; ++it) {
			if (setVolume.find(*it) != setVolume.end()) {
				subVolume.insert(*it);
			}
		}
		scalar++;
	}
	
	return subVolume;
}

void extendSubVolume(std::vector<Z3i::Point>& points,
					 const Z3i::DigitalSet& reference,
					 const Z3i::Point& point,
					 const Z3i::RealPoint& directionVector,
					 double distance) {
	for (auto it = reference.begin(), ite = reference.end();
		 it != ite; ++it) {
		Z3i::Point p(*it);
		if (VCMUtil::abovePlane(p, directionVector, point) && Z3i::l2Metric(p, point) <= distance)
			points.push_back(p);			
	}
}

Z3i::DigitalSet computeNearestPoints(const Z3i::DigitalSet& traversed, const Z3i::Object26_6& reference) {
	Z3i::DigitalSet closestPoints(traversed.domain());
	
	for (auto it = traversed.begin(), ite = traversed.end(); it != ite; ++it) {
		double distanceMin = numeric_limits<double>::max();
		Z3i::Point closest;
		for (auto itO = reference.pointSet().begin(), itOe = reference.pointSet().end(); itO != itOe; ++itO) {
			double distance = Z3i::l2Metric(*itO, *it);
			if (distance < distanceMin) {
				distanceMin = distance;
				closest = *itO;
			}
		}
		closestPoints.insert(closest);
	}
	return closestPoints;
}

template <typename WeightedPointCount>
set<WeightedPointCount*> computeAccumulator(const Z3i::DigitalSet& traversed, const Z3i::DigitalSet& nearest) {
	set<WeightedPointCount*> aSet;
	for (auto it = nearest.begin(), ite = nearest.end(); it != ite; ++it) {
		aSet.insert(new WeightedPointCount(*it, 0, traversed.size()));
	}
	
	for (auto it = traversed.begin(), ite = traversed.end(); it != ite; ++it) {
		for (auto itN = nearest.begin(), itNe = nearest.end(); itN != itNe; ++itN) {
			double distance = Z3i::l2Metric(*it, *itN);
			auto nearestInSetIterator = find_if(aSet.begin(), aSet.end(), [&](WeightedPointCount* wpc) {
					return wpc->myPoint == *itN;
				});
			(*nearestInSetIterator)->myWeight += distance;
		}
	}
	return aSet;
}

template <typename DTL2>
Z3i::DigitalSet computeBallAroundVector(const std::vector<Z3i::Point>& points, const Z3i::DigitalSet& setVolume,
	const DTL2& dt) {
	Z3i::DigitalSet areaOfComputation(setVolume.domain());
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		Ball<Z3i::Point> ball(*it, dt(*it));
		vector<Z3i::Point> pointsInBall = ball.pointsInBall();
		for (auto it = pointsInBall.begin(), ite = pointsInBall.end(); it != ite; ++it) {
			if (setVolume.find(*it) != setVolume.end()) {
				areaOfComputation.insert(*it);
			}
		}
	}
	return areaOfComputation;
}

int minCoordinate(const Z3i::RealPoint& point) {
	int minCoordinate = 0;
	for (unsigned int i = 1; i < Z3i::RealPoint::dimension; i++) {
		if (std::abs(point[i]) < std::abs(point[minCoordinate]))
			minCoordinate = i;
	}
	return minCoordinate;
}

template <typename WeightedPointCount>
WeightedPointCount candidate(const set<WeightedPointCount*>& accumulator, const Z3i::Point& point) {
	double valMax = (*max_element(accumulator.begin(), accumulator.end(), [&](const WeightedPointCount* a, const WeightedPointCount* b) {
				return a->getWeightedValue() < b->getWeightedValue();
			}))->getWeightedValue();

	double valMin = (*min_element(accumulator.begin(), accumulator.end(), [&](const WeightedPointCount* a, const WeightedPointCount* b) {
				return a->getWeightedValue() < b->getWeightedValue();
			}))->getWeightedValue();

	if (valMin == valMax) return WeightedPointCount(point, std::numeric_limits<double>::max(), 1);
	set<WeightedPointCount> reorderedAccumulator;
	for (auto it = accumulator.begin(), ite = accumulator.end(); it != ite; ++it) {
		double value = ((*it)->getWeightedValue() - valMin) / (valMax - valMin);
		double weight = value * Z3i::l2Metric(point, (*it)->myPoint);
		reorderedAccumulator.insert(WeightedPointCount((*it)->myPoint, weight, 1));
	}
	trace.info() << (*(--reorderedAccumulator.end())).myWeight << " " << (*(reorderedAccumulator.begin())).myWeight << endl;
	return (*(--reorderedAccumulator.end()));
}


template <typename DTL2>
Z3i::Point findMaxDTInSet(const Z3i::DigitalSet& set, const DTL2 dt, const Z3i::Point& junctionPoint) {
	double maxDT = 0.0;
	Z3i::Point maxDTPoint;
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		Z3i::Point point = *it;
		if (dt(point) > maxDT) {
			maxDT = dt(point);
			maxDTPoint = point;
		}
		else if (dt(point) == maxDT) {
			if (Z3i::l2Metric(junctionPoint, point) < Z3i::l2Metric(junctionPoint, maxDTPoint))
				maxDTPoint = point;
		}
	}
	return maxDTPoint;
}


vector<Z3i::Point> findEndPoints(const Z3i::DigitalSet& set) {
	Z3i::Object6_26 objectSet(Z3i::dt6_26, set);
	vector<Z3i::Point> endPoints;
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		vector<Z3i::Point> neighbors;
		back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
		objectSet.writeNeighbors(inserter, *it);
		if (neighbors.size() == 1)
			endPoints.push_back(*it);
	}
	return endPoints;
}


template <typename DTL2>
Z3i::Point findLocalMaxDTInSet(const Z3i::DigitalSet& set, const DTL2 dt, const Z3i::Point& junctionPoint) {
	Z3i::Point maxDTPoint;
	double distanceToJunctionPoint = numeric_limits<double>::max();
	Z3i::Object26_6 objSet(Z3i::dt26_6, set);
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		Z3i::Point point = *it;
		
		std::vector<Z3i::Point> neighbors;
		std::back_insert_iterator<std::vector<Z3i::Point>> inserter(neighbors);
		objSet.writeNeighbors(inserter, point);
		
		double distance = dt(point);
		bool candidate = true;
		for (auto itN = neighbors.begin(), itNe = neighbors.end(); itN != itNe; ++itN) {
			if (dt(*itN) >= distance)
				candidate = false;
		}
		if (candidate && Z3i::l2Metric(junctionPoint, point) < distanceToJunctionPoint) {
			maxDTPoint = point;
			distanceToJunctionPoint = Z3i::l2Metric(junctionPoint, point);
		}
	}
	return maxDTPoint;
}

int main( int  argc, char**  argv )
{
	srand(time(NULL));
	QApplication application(argc,argv);
	Viewer3D<> viewer;
	viewer.show();
 
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
	typedef WeightedPoint<Z3i::Point> WeightedPoint;
	typedef WeightedPointCount<Z3i::Point> WeightedPointCount;
	typedef Z3i::KSpace KSpace;
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef DGtal::functors::NotPointPredicate<ThresholdedImage> BackgroundPredicate;
	typedef ImplicitDigitalSurface< KSpace, BackgroundPredicate > DigitalSurfaceContainer;
	
	typedef DGtal::DistanceTransformation<Space, ThresholdedImage, Z3i::L2Metric> DTL2;
	typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric, KernelFunction, DTL2 > VCMOnSurface;
	typedef KSpace::Surfel Surfel;
	typedef Eigen::MatrixXd MatrixXd;
  

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (corresponding volume)")
		("volume,v", po::value<std::string>(), "vol file (corresponding volume)")
		("delta,d", po::value<int>()->default_value(1), "delta")
		("radius,R", po::value<int>()->default_value(5), "big radius")
		("output,o", po::value<std::string>(), "output skeleton filename")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		("thresholdFeature,T",po::value<double>()->default_value(0.1), "threshold feature")
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
	string outFilename = vm["output"].as<std::string>();
	string inputFilename = vm["input"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	int R = vm["radius"].as<int>();
	int delta = vm["delta"].as<int>();
	string volumeFilename = vm["volume"].as<std::string>();
	double thresholdFeature = vm["thresholdFeature"].as<double>();
	Image skeleton = VolReader<Image>::importVol(inputFilename);
	Image volume = VolReader<Image>::importVol(volumeFilename);
	Z3i::DigitalSet setVolume(volume.domain());
	Z3i::Domain domainVolume = setVolume.domain();
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume, 
												  thresholdMin-1, thresholdMax);
	Z3i::Object26_6 object(Z3i::dt26_6, setVolume);
	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);
	set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount> > setVolumeWeighted;
	Image volumeBinary(volume.domain());
	for (auto it = volume.domain().begin(), ite = volume.domain().end(); it != ite; ++it) {
		if (volume(*it) >= thresholdMin && volume(*it) <= thresholdMax) 
			volumeBinary.setValue(*it, 255);
		else
			volumeBinary.setValue(*it, 0);
	}

	Z3i::DigitalSet skeletonPoints(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image>(skeletonPoints, skeleton,
												 thresholdMin-1, thresholdMax);
		
   

	ThresholdedImage binarizer(volume, thresholdMin-1, thresholdMax);
	DTL2 dt(&volume.domain(), &binarizer, &Z3i::l2Metric);
	BackgroundPredicate backgroundPredicate(binarizer);

	DTL2::Domain domainDT = dt.domain();
	for (auto it = domainDT.begin(), ite = domainDT.end(); it != ite; ++it) {
		double value = dt(*it);
		if (value > 0) {
			setVolumeWeighted.insert(new WeightedPointCount(*it, value));
		}		
	}

	double distanceMax = (*setVolumeWeighted.begin())->myWeight;
	KSpace ks;
	Metric l2;
	
	ks.init( volume.domain().lowerBound(),
			 volume.domain().upperBound(), true );
	SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
	Surfel bel = Surfaces<KSpace>::findABel( ks, backgroundPredicate, 1000000 );
	DigitalSurfaceContainer* container =
		new DigitalSurfaceContainer( ks, backgroundPredicate, surfAdj, bel, false  );
	DigitalSurface< DigitalSurfaceContainer > surface( container ); //acquired
	
	//! [DVCM3D-instantiation]
	Surfel2PointEmbedding embType = Pointels; // Could be Pointels|InnerSpel|OuterSpel;
	KernelFunction chiSurface( 1.0, R );             // hat function with support of radius r
	VCMOnSurface* vcm_surface = new VCMOnSurface( surface, embType, R, 1, chiSurface, dt, delta, l2, true);
	Z3i::DigitalSet branchingPoints = computeBranchingPartsWithVCMFeature(*vcm_surface, domainVolume, thresholdFeature);
	Z3i::Object26_6 obj(Z3i::dt26_6, branchingPoints);
	vector<Z3i::Object26_6> objectsBranching;
	back_insert_iterator< std::vector<Z3i::Object26_6> > inserterBranching( objectsBranching );
	obj.writeComponents(inserterBranching);
	Z3i::DigitalSet maxCurvaturePoints(domainVolume);
	for (auto it = objectsBranching.begin(), ite = objectsBranching.end(); it != ite; ++it) {
		double ratioMax = 0;
		Point maximizingCurvaturePoint;
		for (auto itPoint = it->pointSet().begin(), itPointE = it->pointSet().end(); itPoint != itPointE; ++itPoint) {
			auto lambda = (vcm_surface->mapPoint2ChiVCM()).at(*itPoint).values;
			double ratio = lambda[0] / (lambda[0] + lambda[1] + lambda[2]); 
			if (ratio > ratioMax) {
				ratioMax = ratio;
				maximizingCurvaturePoint = *itPoint;
			}
		}
   
 		maxCurvaturePoints.insert(maximizingCurvaturePoint);
	}
	
	

	
	VCM vcm( R, 1, l2, true );
	vcm.init( setVolume.begin(), setVolume.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, 1 );
	
	typedef Z3i::Object26_6 ObjectType;


	set<Z3i::Point> pointsToProcess;

	ObjectType objectImage(Z3i::dt26_6, skeletonPoints);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);
	trace.info() << nbConnectedComponents << endl;
	sort(objects.begin(), objects.end(), [&](ObjectType a, ObjectType b)->bool {
			double maxa = 0, maxb = 0;
			for (auto it = a.begin(), ite = a.end(); it != ite; ++it) {
				if (dt(*it) > maxa)
					maxa = dt(*it);
			}
			
			for (auto it = b.begin(), ite = b.end(); it != ite; ++it) {
				if (dt(*it) > maxb)
					maxb = dt(*it);
			}
			return maxa > maxb;
		});
	
	ObjectType reference = *(objects.begin());
	objects.erase(objects.begin());

	Point closestBranchingCurrent;
	trace.beginBlock("Connecting disconnected components");
   
	while (objects.size() > 0) {
		trace.progressBar(nbConnectedComponents-objects.size(), nbConnectedComponents);
		vector<ObjectType>::iterator minimizingObjectToReference = objects.end();
		double distanceMin = numeric_limits<double>::max();
		vector<Z3i::Point> points;
		Z3i::Point belongingToCurrentObject, belongingToReference;
//		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
		ObjectType currentObject = *(objects.begin());
		
		for (auto itCurrentObject = currentObject.pointSet().begin(),
				 itCurrentObjectE = currentObject.pointSet().end(); itCurrentObject != itCurrentObjectE; ++itCurrentObject) {

			for (auto itReference = reference.pointSet().begin(), itReferenceE = reference.pointSet().end(); itReference != itReferenceE; ++itReference) {
				Point closestBranchingPoint = *min_element(maxCurvaturePoints.begin(), maxCurvaturePoints.end(), [&](const Point& one, const Point& two) {
						return Z3i::l2Metric(one, *itCurrentObject) < Z3i::l2Metric(two, *itCurrentObject);
					});
				vector<Point> points = PointUtil::linkTwoPoints(*itReference, *itCurrentObject);
				bool add = true;
				double distanceDT = min(dt(*itCurrentObject), dt(*itReference));
				for (auto itP = points.begin(), itPe = points.end(); itP != itPe; ++itP) {
					if (Z3i::l2Metric(*itP, closestBranchingPoint)+1 < distanceDT ||
						dt(*itP) == 0) {
						add = false;
						break;
					}
				}
					
				double distance = Z3i::l2Metric(*itCurrentObject, *itReference);
				if (distance < distanceMin && add) {
					minimizingObjectToReference = objects.begin();
					distanceMin = distance;
					belongingToCurrentObject = *itCurrentObject;
					belongingToReference = *itReference;
					closestBranchingCurrent = closestBranchingPoint;
				}
			}
		}
		//}
		if (minimizingObjectToReference == objects.end()) {
			minimizingObjectToReference = objects.begin();
			objects.erase(minimizingObjectToReference);
		}
		else if (Z3i::l2Metric(belongingToReference, belongingToCurrentObject) <= 2 * sqrt(3)) {
			for (auto it = minimizingObjectToReference->pointSet().begin(), ite = minimizingObjectToReference->pointSet().end();
				 it != ite; ++it) {
				reference.pointSet().insert(*it);
			}
			vector<Point> newPoints = PointUtil::linkTwoPoints(belongingToCurrentObject, belongingToReference);
			for (auto it = newPoints.begin(), ite = newPoints.end(); it != ite; ++it) 
				skeletonPoints.insert(*it);
			objects.erase(minimizingObjectToReference);
		}
		else {
	  
		   
			// vector<Point> newPoints = PointUtil::linkTwoPoints(belongingToCurrentObject, belongingToReference);
			// Z3i::RealPoint normal = (belongingToReference - belongingToCurrentObject).getNormalized();
			// Point current = belongingToCurrentObject;
			// Z3i::DigitalSet curve(domain);
			// do { 
			// 	vector<Point> neighbors;
			// 	back_insert_iterator<vector<Point>> inserter(neighbors);
			// 	MetricAdjacency<Space, 3>::writeNeighbors(inserter, current);
			// 	Z3i::DigitalSet setNeighbors(domain);
			// 	for (auto it = neighbors.begin(), ite = neighbors.end(); it != ite; ++it) {
			// 		if (VCMUtil::abovePlane(*it, normal, current) && curve.find(*it) == curve.end())
			// 			setNeighbors.insert(*it);
			// 	}
			// 	Z3i::Point p = findMaxDTInSet(setNeighbors, dt, current);
			// 	viewer << CustomColors3D(Color::Blue, Color::Blue) << p;
			// 	current = p;
			// 	curve.insert(current);
			// } while (VCMUtil::abovePlane(belongingToReference, normal, current));
			
			// viewer << Viewer3D<>::updateDisplay;
			// qApp->processEvents();
			for (auto it = minimizingObjectToReference->pointSet().begin(), ite = minimizingObjectToReference->pointSet().end();
 				 it != ite; ++it) {
 				reference.pointSet().insert(*it);
 			}
			
			double radiusCurrentObject = dt(belongingToCurrentObject) ;
			chi = KernelFunction( 1.0, radiusCurrentObject );
			vcm.setMySmallR(radiusCurrentObject);
			Z3i::RealPoint normalCurrentObject;
			VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
										  belongingToCurrentObject, normalCurrentObject,
										  0, radiusCurrentObject, distanceMax, true);
			

			double radiusReference = dt(belongingToReference);
			chi = KernelFunction( 1.0, radiusReference );
			vcm.setMySmallR(radiusReference);
			Z3i::RealPoint normalReference;
			VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
										  belongingToReference, normalReference,
										  0, radiusReference, distanceMax, true);
			Z3i::RealPoint dirVectorCurrent = (belongingToReference - belongingToCurrentObject).getNormalized();
			if (normalCurrentObject.dot(dirVectorCurrent) < 0)
				normalCurrentObject = -normalCurrentObject;

			Z3i::RealPoint dirVectorReference = (belongingToCurrentObject - belongingToReference).getNormalized();
			if (normalReference.dot(dirVectorReference) < 0)
				normalReference = -normalReference;
			
			vector<Point> newPoints;

			
//			if (Z3i::l2Metric(belongingToReference, belongingToCurrentObject) <= 2 * sqrt(3))
			newPoints = PointUtil::linkTwoPoints(belongingToCurrentObject, belongingToReference);
//			else	
//				newPoints = PointUtil::bezierCurve(belongingToCurrentObject, belongingToReference, controlCurrent, controlReference);

			vector<Point> computationPoints = newPoints;
			extendSubVolume(newPoints, reference.pointSet(), belongingToReference, -normalReference, dt(belongingToReference));		   
			extendSubVolume(newPoints, reference.pointSet(), belongingToCurrentObject, -normalCurrentObject, dt(belongingToCurrentObject));
		 	  
			Z3i::DigitalSet newPointsSet(setVolume.domain());
			newPointsSet.insert(newPoints.begin(), newPoints.end());
			for (auto it = newPoints.begin(), ite = newPoints.end(); it != ite; ++it) {
				//skeletonPoints.insert(*it);
				//viewer << CustomColors3D(Color::Blue, Color::Blue) << *it;
			}
			// viewer << Viewer3D<>::updateDisplay;
			// qApp->processEvents();

			Z3i::DigitalSet computationVolume = computeBallAroundVector(newPoints, setVolume, dt);
			Z3i::DigitalSet restrictedComputationVolume = computeBallAroundVector(computationPoints, setVolume, dt);
			int r = rand() % 256, g = rand() % 256, b = rand() % 256;
			Color c(r,g,b);
//			viewer << CustomColors3D(c,c) << computationVolume;
			
			set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount> > subVolumeWeighted;
			for (auto it = computationVolume.begin(), ite = computationVolume.end(); it != ite; ++it) {
				auto itW = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
						return (wpc->myPoint == *it);
					});
				if (itW != setVolumeWeighted.end()) {
					subVolumeWeighted.insert(new WeightedPointCount(*(*itW)));
				}
			}

			set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount> > subVolumeWeightedComputation;
			for (auto it = restrictedComputationVolume.begin(), ite = restrictedComputationVolume.end(); it != ite; ++it) {
				auto itW = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
						return (wpc->myPoint == *it);
					});
				if (itW != setVolumeWeighted.end()) {
					subVolumeWeightedComputation.insert(new WeightedPointCount(*(*itW)));
				}
			}
			
			Z3i::DigitalSet connectedComponent3D(setVolume.domain());
			Z3i::RealPoint realCenter, normalSub;
			
			WeightedPointCount* currentWPC;
			double distance = std::numeric_limits<double>::max();
			for (auto it = subVolumeWeightedComputation.begin(), ite = subVolumeWeightedComputation.end(); it != ite; ++it) {
				double currentDistance = Z3i::l2Metric((*it)->myPoint, belongingToCurrentObject);				
				if (currentDistance < distance) {
					distance = currentDistance;
					currentWPC = *it;
				}
			}
			int numberLeft = subVolumeWeightedComputation.size();
			while (numberLeft > 0 && (VCMUtil::abovePlane(belongingToReference, normalSub, realCenter))) {
//			for (auto it = newPointsSet.begin(), ite = newPointsSet.end(); it != ite; ++it) {
				// auto currentPoint = find_if(subVolumeWeighted.begin(), subVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
				// 		return wpc->myPoint == *it;
				// 	});
				// if (currentPoint == subVolumeWeighted.end()) continue;
				// currentWPC = *currentPoint;
				
				currentWPC->myProcessed = true;
				double radius = dt(currentWPC->myPoint);
				connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, subVolumeWeighted, currentWPC->myPoint, normalSub, 0,  radius, 100, true);
				
				realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);
				Z3i::Point centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
				viewer << CustomColors3D(Color::Blue, Color::Blue) << centerOfMass;
				viewer << Viewer3D<>::updateDisplay;
				qApp->processEvents();
				bool processed = false;
				for (auto it = connectedComponent3D.begin(), ite = connectedComponent3D.end(); it != ite; ++it) {
					for (auto itS = skeletonPoints.begin(), itSe = skeletonPoints.end(); itS != itSe; ++itS)  {
						if (*itS == *it)
							processed = true;
					}
				}
				if (!processed &&  Z3i::l2Metric(currentWPC->myPoint, centerOfMass) <= sqrt(3))
					skeletonPoints.insert(centerOfMass);

				if (normalSub.dot(dirVectorCurrent) < 0)
					normalSub = -normalSub;
				VCMUtil::markConnectedComponent3D(subVolumeWeightedComputation, connectedComponent3D, 0);
				VCMUtil::trackNextPoint(currentWPC, subVolumeWeightedComputation, connectedComponent3D, centerOfMass, normalSub);
				numberLeft = count_if(subVolumeWeightedComputation.begin(), subVolumeWeightedComputation.end(),
									  [&](WeightedPointCount* wpc) {
										  return (!wpc->myProcessed);
									  });
				if (currentWPC == (*subVolumeWeightedComputation.end()))
					break;
			}
			
			objects.erase(minimizingObjectToReference);			
		}
	}
	trace.endBlock();

	delete vcm_surface;
	for (auto it = skeletonPoints.begin(), ite = skeletonPoints.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color::Red, Color::Red) << *it;
	}
	for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color(0,0,255,30), Color(0,0,255,30)) << *it;
	}
	
	Image outImage(volume.domain());
	
	DGtal::imageFromRangeAndValue(skeletonPoints.begin(), skeletonPoints.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);

	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
	   
