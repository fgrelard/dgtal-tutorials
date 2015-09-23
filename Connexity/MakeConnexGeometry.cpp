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

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


int main( int  argc, char**  argv )
{
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

Image skeleton = VolReader<Image>::importVol(inputFilename);
Image volume = VolReader<Image>::importVol(volumeFilename);
Z3i::DigitalSet setVolume(volume.domain());
Z3i::Domain domainVolume = setVolume.domain();
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume, 
												  thresholdMin-1, thresholdMax);
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
	
  
   
	typedef DGtal::functors::IntervalForegroundPredicate<Image> Binarizer; 
	Binarizer binarizer(volume, thresholdMin-1, thresholdMax);
	typedef DGtal::DistanceTransformation<Space, Binarizer, Z3i::L2Metric> DTL2;
	
	DTL2 dt(&volume.domain(), &binarizer, &Z3i::l2Metric);
		

	DTL2::Domain domainDT = dt.domain();
	for (auto it = domainDT.begin(), ite = domainDT.end(); it != ite; ++it) {
		double value = dt(*it);
		if (value > 0) {
			setVolumeWeighted.insert(new WeightedPointCount(*it, value));
		}		
	}

	for (auto it = skeletonPoints.begin(), ite = skeletonPoints.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color::Red, Color::Red) << *it;
	}
	
	Metric l2;
	VCM vcm( R, R, l2, true );
	vcm.init( setVolume.begin(), setVolume.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, R );
	
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
	trace.beginBlock("Connecting disconnected components");
	while (objects.size() > 0) {
		trace.progressBar(nbConnectedComponents-objects.size(), nbConnectedComponents);
		vector<ObjectType>::iterator minimizingObjectToReference;
		double distanceMin = numeric_limits<int>::max();
		Z3i::Point belongingToCurrentObject, belongingToReference;
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			ObjectType currentObject = *it;
			for (auto itCurrentObject = currentObject.pointSet().begin(),
					 itCurrentObjectE = currentObject.pointSet().end(); itCurrentObject != itCurrentObjectE; ++itCurrentObject) {				
				for (auto itReference = reference.pointSet().begin(), itReferenceE = reference.pointSet().end(); itReference != itReferenceE; ++itReference) {
					Z3i::RealPoint vectorDifference = *itReference - *itCurrentObject;
					double distance = Z3i::l2Metric(*itCurrentObject, *itReference);
					if (distance < distanceMin) {
						minimizingObjectToReference = it;
						distanceMin = distance;
						belongingToCurrentObject = *itCurrentObject;
						belongingToReference = *itReference;
					}
				}
			}
		}

		double radius = dt(belongingToCurrentObject) + delta;
		vcm.updateProximityStructure(radius, setVolume.begin(), setVolume.end());
		chi = KernelFunction( 1.0, radius);

		Z3i::RealPoint normal = VCMUtil::computeNormalFromVCM(belongingToCurrentObject, vcm, chi, 0);

		Z3i::Point projection(belongingToCurrentObject);
		Z3i::Point previousProjection(belongingToCurrentObject);
		int scalar = 1;
	    vector<WeightedPointCount> accumulator;
		while (reference.pointSet().find(projection) == reference.pointSet().end() &&
			   setVolume.find(projection) != setVolume.end()) {

			projection = belongingToCurrentObject + normal * scalar;
			
			if (projection != previousProjection) {
				double distanceMin = numeric_limits<int>::max();
				Z3i::Point referencePoint;
				for (auto it = reference.pointSet().begin(), ite = reference.pointSet().end();
					 it != ite; ++it) {
					double distance = Z3i::l2Metric(*it, projection);
					if (distance < distanceMin) {
						distanceMin = distance;
						referencePoint = *it;
					}
				}
				auto itAccumulator = find_if(accumulator.begin(), accumulator.end(), [&](const WeightedPointCount& wpc) {
						return wpc.myPoint == referencePoint;
					});
				if (itAccumulator  == accumulator.end()) {
					accumulator.push_back(WeightedPointCount(referencePoint, distanceMin, 1));
				} else {
					itAccumulator->myWeight += distanceMin;
					itAccumulator->myCount++;
				}
			}
		
			scalar++;	
			previousProjection = projection;
			
		}

		sort(accumulator.begin(), accumulator.end(), [](const WeightedPointCount& a, const WeightedPointCount& b)
			 {
				 return a.getWeightedValue() < b.getWeightedValue();
			 }
			);
		WeightedPointCount propositionOne( *(accumulator.begin()));
		
		accumulator.clear();
		projection = belongingToCurrentObject;
		scalar = 1;
		while (reference.pointSet().find(projection) == reference.pointSet().end() &&
			   setVolume.find(projection) != setVolume.end()) {
			projection = belongingToCurrentObject - normal * scalar;
			
			if (projection != previousProjection) {
				double distanceMin = numeric_limits<int>::max();
				Z3i::Point referencePoint;
				for (auto it = reference.pointSet().begin(), ite = reference.pointSet().end();
					 it != ite; ++it) {
					double distance = Z3i::l2Metric(*it, projection);
					if (distance < distanceMin) {
						distanceMin = distance;
						referencePoint = *it;
					}
				}
				auto itAccumulator = find_if(accumulator.begin(), accumulator.end(), [&](const WeightedPointCount& wpc) {
						return wpc.myPoint == referencePoint;
					});
				if (itAccumulator  == accumulator.end()) {
					accumulator.push_back(WeightedPointCount(referencePoint, distanceMin, 1));
				} else {
					itAccumulator->myWeight += distanceMin;
					itAccumulator->myCount++;
				}
			}
		
			scalar++;	
			previousProjection = projection;
		}
		
		sort(accumulator.begin(), accumulator.end(), [](const WeightedPointCount& a, const WeightedPointCount& b)
			 {
				 return a.getWeightedValue() < b.getWeightedValue();
			 }
			);
	
		WeightedPointCount propositionTwo(*(accumulator.begin()));
		
		if (propositionOne.getWeightedValue() < propositionTwo.getWeightedValue()) {
			belongingToReference = propositionOne.myPoint;
		} else  {
			belongingToReference = propositionTwo.myPoint;
		}
		for (auto it = minimizingObjectToReference->pointSet().begin(),
				 ite = minimizingObjectToReference->pointSet().end();
			 it != ite; ++it) {
			reference.pointSet().insert(*it);
		}
		
		
		vector<Z3i::Point> points = PointUtil::linkTwoPoints(belongingToCurrentObject, belongingToReference);
		for (auto it = points.begin(), ite = points.end(); it!=ite; ++it) {
			viewer << CustomColors3D(Color::Blue, Color::Blue) << *it;
			viewer << Viewer3D<>::updateDisplay;
			qApp->processEvents();
			skeletonPoints.insert(*it);
			pointsToProcess.insert(*it);
			reference.pointSet().insert(*it);
		}
		objects.erase(minimizingObjectToReference);
	}

	
	
	// double radius = 0;
	// Z3i::RealPoint normal, previousNormal;
	// int i = 0;
	// for (auto it = pointsToProcess.begin(), ite = pointsToProcess.end();
	// 	 it != ite; ++it) {
	// 	trace.progressBar(i, pointsToProcess.size());
	// 	radius = dt(*it) + delta;
	// 	vcm.updateProximityStructure(radius, setVolume.begin(), setVolume.end());
	// 	chi = KernelFunction( 1.0, radius);
	// 	Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
	// 																*it, normal, 0, previousNormal, radius);
	// 	Z3i::Point centerOfMass = Statistics::extractCenterOfMass3D(connectedComponent3D);
				
	// 	if (centerOfMass != Z3i::Point()) {
	// 		reference.pointSet().insert(centerOfMass);
	// 		skeletonPoints.insert(centerOfMass);
	// 		viewer << CustomColors3D(Color::Blue, Color::Blue) << centerOfMass;
	// 		viewer << Viewer3D<>::updateDisplay;
	// 		qApp->processEvents();
	// 	}
	// 	i++;
	// }
	trace.endBlock();
	


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
	   
