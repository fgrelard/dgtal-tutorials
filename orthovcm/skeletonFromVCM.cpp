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
#include <queue>
#include <iterator>
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
#include "WeightedPointCount.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

class LabelledPoint : public WeightedPoint<Z3i::Point> {
	typedef WeightedPoint<Z3i::Point> Base;
	using Base::Base;
	friend bool operator<(const LabelledPoint& it, const LabelledPoint& other) {
		return (it.myPoint < other.myPoint);
	}
};

struct WeightedPointCountComparator
{
	bool operator()(const WeightedPointCount<Z3i::Point>* lhs, const WeightedPointCount<Z3i::Point>* rhs) const  {
		return lhs->myWeight >= rhs->myWeight;
	}
};


template <typename Image>
Z2i::DigitalSet extractConnectedComponent(const Image& image, const Z2i::Point& referencePoint, int thresholdMin,
										  int thresholdMax) {
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image2D;
	typedef Z2i::Object8_4 ObjectType;
	
	//Extracting connected components
	Z2i::DigitalSet points2D(image.domain());
	SetFromImage<Z2i::DigitalSet>::append<Image> (points2D, image, 
												  thresholdMin-1, thresholdMax);
	ObjectType objectImage2D(Z2i::dt8_4, points2D);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector< ObjectType > > inserter( objects );
	unsigned int nbConnectedComponents = objectImage2D.writeComponents(inserter);
	Image2D image2D(image.domain());
	if (nbConnectedComponents > 1) {
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			Z2i::DigitalSet ccSet = it->pointSet();
			if (ccSet.find(referencePoint) != ccSet.end()) {
			    return ccSet;
			}
		}
		return Z2i::DigitalSet(image.domain());
	}
	return points2D;
}

template <typename Domain, typename WeightedPoint>
Z3i::DigitalSet extractConnectedComponent3D(const Domain & domain, const set<WeightedPoint*, WeightedPointCountComparator>& volume, const Z3i::RealPoint& normal, const Z3i::Point& referencePoint, double d, double omega, double distanceMax) {
	typedef Z3i::Object26_6 ObjectType;

	Z3i::DigitalSet intersection(domain);
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		double valueToCheckForPlane = (*it)->myPoint[0] * normal[0] + (*it)->myPoint[1] * normal[1] + (*it)->myPoint[2] * normal[2];
		if (valueToCheckForPlane >= d && valueToCheckForPlane < d + omega && Z3i::l2Metric((*it)->myPoint, referencePoint) <= distanceMax) {
			intersection.insert((*it)->myPoint);
		}
	}

	ObjectType objectIntersection(Z3i::dt26_6, intersection);
	vector<ObjectType> objects;
	back_insert_iterator<vector<ObjectType>> inserter(objects);
	unsigned int nbConnectedComponents = objectIntersection.writeComponents(inserter);
	if (nbConnectedComponents > 1) {
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			Z3i::DigitalSet ccSet = it->pointSet();
			if (ccSet.find(referencePoint) != ccSet.end()) {
			    return ccSet;
			}
		}
		return Z3i::DigitalSet(domain);
	}
	return intersection;
}

template <typename WeightedPoint>
bool markConnectedComponent3D(set<WeightedPoint*, WeightedPointCountComparator>& volume, const Z3i::DigitalSet& intersection, int label) {	

	set<int> processedLabels;
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		if (!(*it)->myProcessed && intersection.find((*it)->myPoint) != intersection.end()) {
			(*it)->myProcessed = true;
			(*it)->myCount = label;
		}
		else if ((*it)->myProcessed  && intersection.find((*it)->myPoint) != intersection.end()) {
			processedLabels.insert((*it)->myCount);
		}
	}
	return (processedLabels.size() <= 3);	
}

template <typename Adapter>
Z3i::DigitalSet project2DSetIn3D(const Z2i::DigitalSet& set, const Z3i::Domain& domain3D, const Adapter& adapter) {
	Z3i::DigitalSet set3D(domain3D);
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		set3D.insert(adapter(*it));
	}
	return set3D;
}

Z2i::RealPoint extractCenterOfMass(const Z2i::DigitalSet& set) {
	if (set.size() != 0) {
		typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image2D;
		Image2D image2D = ImageFromSet<Image2D>::create(set, 150);
		Z2i::RealPoint centerOfMass = SliceUtils::centerOfMass(image2D);
		return centerOfMass;
	}
	else
		return Z2i::RealPoint();
}

Z3i::RealPoint extractCenterOfMass3D(const Z3i::DigitalSet& set) {
	if (set.size() != 0) {
		typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image3D;
		Image3D image3D = ImageFromSet<Image3D>::create(set, 150);
		Z3i::RealPoint centerOfMass = SliceUtils::centerOfMass3D(image3D);
		return centerOfMass;
	}
	else
		return Z3i::RealPoint();
}

Z3i::DigitalSet computePointsBelowPlane(const Z3i::DigitalSet& volume, const Z3i::RealPoint& normal, const Z3i::Point& center) {
	Z3i::DigitalSet inferiorPoints(volume.domain());
	double d = -(-normal[0] * center[0] - normal[1] * center[1] - normal[2] * center[2]);
	for (auto it = volume.begin(), ite = volume.end(); it!=ite; ++it) {
		double valueToCheckForPlane = (*it)[0] * normal[0] + (*it)[1] * normal[1] + (*it)[2] * normal[2];
		if (valueToCheckForPlane < d)
			inferiorPoints.insert(*it);
	}
	return inferiorPoints;
}

set<Z3i::Point> differenceBetweenPlanes(const Z3i::DigitalSet& volume, const Z3i::RealPoint& normalPrevious, const Z3i::Point& centerPrevious,
										const Z3i::RealPoint& normalNext, const Z3i::Point& centerNext) {
	Z3i::RealPoint next = normalNext;
	if (normalPrevious.dot(normalNext) < 0) {
		next  = -normalNext;
 	}
	Z3i::DigitalSet inferiorPointsPrevious = computePointsBelowPlane(volume, normalPrevious, centerPrevious);
	Z3i::DigitalSet inferiorPointsNext = computePointsBelowPlane(volume, next, centerNext);
		
	set<Z3i::Point> difference;
	if (inferiorPointsNext.size() > inferiorPointsPrevious.size()) {
		set_difference(inferiorPointsNext.begin(), inferiorPointsNext.end(),
					   inferiorPointsPrevious.begin(), inferiorPointsPrevious.end(),
					   inserter(difference, difference.end()));
	}
	else {
		set_difference(inferiorPointsPrevious.begin(), inferiorPointsPrevious.end(),
					   inferiorPointsNext.begin(), inferiorPointsNext.end(),
					   inserter(difference, difference.end()));
	}
	return difference;
}

template <typename WeightedPoint>
void computeDistanceFromCenterOfMass(set<WeightedPoint*>& weightedSet, const Z3i::RealPoint& centerOfGravity) {
	
	Z3i::Point discreteCoG(round(centerOfGravity[0]), round(centerOfGravity[1]), round(centerOfGravity[2]));
	
	vector<Z3i::Point> neighbors;
	back_insert_iterator<vector<Z3i::Point>> iterator(neighbors);
	MetricAdjacency<Z3i::Space, 3>::writeNeighbors(iterator, discreteCoG);

	double distance_max = 0;
	double distance_min = numeric_limits<double>::max();
	vector<WeightedPoint*> neighborWeighted;
	for (auto& it : weightedSet) {
		if (find(neighbors.begin(), neighbors.end(), it->myPoint) != neighbors.end()) {
			double d = sqrt(pow(centerOfGravity[0] - it->myPoint[0], 2) +
							pow(centerOfGravity[1] - it->myPoint[1], 2) +
							pow(centerOfGravity[2] - it->myPoint[2], 2));
			if (d > distance_max) distance_max = d;
			if (d < distance_min) distance_min = d;
			neighborWeighted.push_back(it);
		}
	}
	

	for (auto& it : neighborWeighted) {
		it->myCount++;
		double d = sqrt(pow(centerOfGravity[0] - it->myPoint[0], 2) +
						pow(centerOfGravity[1] - it->myPoint[1], 2) +
						pow(centerOfGravity[2] - it->myPoint[2], 2));
		it->myWeight += 1- ((d - distance_min) / (distance_max - distance_min));
	}
}

Z3i::DigitalSet computeSetOfSimplePoints(const Z3i::DigitalSet& points) {
	typedef Z3i::Object26_6 ObjectType;

	Z3i::DigitalSet simplePoints(points.domain());
	ObjectType object(Z3i::dt26_6, points);

	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		if (object.isSimple(*it)) {
			simplePoints.insert(*it);
		}
	}

	return simplePoints;
}

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
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef VoronoiCovarianceMeasure<Space,Metric> VCM;
	typedef functors::HatPointFunction<Point,double> KernelFunction;
  
	typedef MSTTangent<Point> Tangent;
	typedef Pencil<Point, Tangent, RealPoint> Pencil;
	typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;
	typedef WeightedPoint<Z3i::RealPoint> WeightedRealPoint;
	typedef WeightedPoint<Z3i::Point> WeightedPoint;
	typedef MetricAdjacency<Space, 3> MetricAdjacency;
	typedef WeightedPointCount<Point> WeightedPointCount;
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;
	
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("skeleton,s", po::value<std::string>(), "vol file (medial axis)")
		("input,i", po::value<std::string>(), "vol file (corresponding volume)")
		("output,o", po::value<std::string>(), "output skeleton filename")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		("radiusInside,R", po::value<double>()->default_value(10), "radius of the ball inside voronoi cell")
		("radiusNeighbour,r", po::value<double>()->default_value(10), "radius of the ball for the neighbourhood")
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
	double R = vm["radiusInside"].as<double>();
	double r = vm["radiusNeighbour"].as<double>();


	Image volume = VolReader<Image>::importVol(inputFilename);
	Z3i::Domain domainVolume = volume.domain();
	Z3i::DigitalSet setVolume(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume, 
												  thresholdMin-1, thresholdMax);
	set<WeightedPointCount*, WeightedPointCountComparator> setVolumeWeighted;
	
	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);

	vector<Point> vPoints;
	vector<Pencil> tangents;
	if (vm.count("skeleton")) {
		string skeletonFilename = vm["skeleton"].as<std::string>();
		Image image = VolReader<Image>::importVol(skeletonFilename);
		Z3i::Domain domainSkeleton = image.domain();
		Z3i::DigitalSet setSkeleton(domainSkeleton);
		SetFromImage<Z3i::DigitalSet>::append<Image> (setSkeleton, image, 
													  thresholdMin-1, thresholdMax);
	
	
		Point p;	
		for (auto it = image.domain().begin(), itE = image.domain().end(); it != itE; ++it) {
			if (image(*it) >= thresholdMin && image(*it) <= thresholdMax) {
				p = *it;
			}
		}
	

		vPoints = PointUtil::containerFromDepthTraversal<vector<Point>>(image, p, thresholdMin, thresholdMax);

		//Computing lambda MST tangents
		tangents = TangentUtils::orthogonalPlanesWithTangents<Pencil>(vPoints.begin(), vPoints.end());
	}
	

	Z3i::DigitalSet skeletonPoints(domainVolume);
	
	typedef DGtal::functors::IntervalForegroundPredicate<Image> Binarizer; 
	Binarizer binarizer(volume, thresholdMin-1, thresholdMax);
	typedef DGtal::DistanceTransformation<Space, Binarizer, Z3i::L2Metric> DTL2;
	
	DTL2 dt(&volume.domain(), &binarizer, &Z3i::l2Metric);

	for (auto it = dt.domain().begin(), ite = dt.domain().end(); it != ite; ++it) {
		double value = dt(*it);
		if (value > 0)
			setVolumeWeighted.insert(new WeightedPointCount(*it, value));
	}
	
	WeightedPointCount* currentPoint = *setVolumeWeighted.begin();
	const double size = 20.0; // size of displayed normals

  
	Metric l2;
	VCM vcm( R, ceil( r ), l2, true );
	vcm.init( setVolume.begin(), setVolume.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, r );


	Matrix vcm_r, evec;
	RealVector eval;
	Viewer3D<> viewer;
	viewer.show();
 
	const int IMAGE_PATCH_WIDTH = 100;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH), volume.domain().upperBound() + Z3i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
									  DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::functors::Identity idV;
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );

	int i = 0;
	int numberLeft = setVolumeWeighted.size();
	Z3i::RealPoint previousNormal;
	Z3i::Point previousCenter;
	Z3i::DigitalSet previousConnectedComponent3D(domainVolume);
	
	trace.beginBlock("Computing skeleton");
	while (numberLeft > 0)
	{
		currentPoint->myProcessed = true;
		trace.progressBar((setVolume.size() - numberLeft), setVolume.size());		
		double radius = r;
		if (vm.count("skeleton")) {
			Pencil closestPointToCurrent = *min_element(tangents.begin(), tangents.end(), [&](const Pencil& one, const Pencil& two) {
					return Z3i::l2Metric(one.getPoint(), currentPoint->myPoint) < Z3i::l2Metric(two.getPoint(), currentPoint->myPoint);
				});
			
			DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain> embedder(domain3Dyup, closestPointToCurrent.getPoint(), closestPointToCurrent.getTangent(), IMAGE_PATCH_WIDTH);
			ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
			extractedImage.setDefaultValue(0);
			radius  = SliceUtils::computeRadiusFromImage(extractedImage, thresholdMin, thresholdMax);
			radius += 2;
			if (radius > 0) {
				vcm.updateProximityStructure(radius*2, setVolume.begin(), setVolume.end());
				chi = KernelFunction( 1.0, radius*2);
			}
		}
		
		// Compute VCM and diagonalize it.
		vcm_r = vcm.measure( chi, currentPoint->myPoint );
		LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );

		// // Display normal
		RealVector normal = evec.column(0);
		if (normal.dot(previousNormal) < 0)
			 normal = -normal;
		RealVector n = evec.column( 2 );
		n*=size;
		RealPoint p( currentPoint->myPoint[ 0 ], currentPoint->myPoint[ 1 ], currentPoint->myPoint[2] );
		RealVector n2 = evec.column( 1 );
		n2*=size;

		viewer.setLineColor(Color::Red);
		viewer.setFillColor(Color::Red);
		viewer.setFillTransparency(150);
		
		double imageSize = radius*2;

		double miniCoordinatesInNormal = min(abs(normal[0]), min(abs(normal[1]), abs(normal[2])));
		double secondMini = min(abs(normal[0]), abs(normal[1])) > miniCoordinatesInNormal ? min(abs(normal[0]), abs(normal[1])) :
			min(abs(normal[0]), abs(normal[2])) > miniCoordinatesInNormal ? min(abs(normal[0]), abs(normal[2])) :
			min(abs(normal[1]), abs(normal[2]));
		double factor = 1. / min((secondMini - miniCoordinatesInNormal), miniCoordinatesInNormal);
		if (factor > 1000)
			factor = 1;
		Z3i::Point discreteNormal(round(normal[0]*factor), round(normal[1]*factor), round(normal[2]*factor));
		double d = -(-discreteNormal[0] * currentPoint->myPoint[0] - discreteNormal[1] * currentPoint->myPoint[1] - discreteNormal[2] * currentPoint->myPoint[2]);
		double omega = abs(discreteNormal[0]) + abs(discreteNormal[1]) + abs(discreteNormal[2]);
	
		Z3i::DigitalSet connectedComponent3D = extractConnectedComponent3D(domainVolume, setVolumeWeighted, discreteNormal, currentPoint->myPoint, d, omega, imageSize);
		Z3i::RealPoint centerOfMass = extractCenterOfMass3D(connectedComponent3D);

//		 for (auto it = connectedComponent3D.begin(), ite = connectedComponent3D.end();
//			  it != ite; ++it) {
//			 viewer << CustomColors3D(Color::Green, Color::Green) << *it;
//		 }
		
		//computeDistanceFromCenterOfMass(setVolumeWeighted, centerOfMassEmbedded);
		
		vector<Point> neighbors26;
		back_insert_iterator<vector<Point>> iterator(neighbors26);
		MetricAdjacency::writeNeighbors(iterator, centerOfMass);
		int cpt = 0;
		for (auto it = neighbors26.begin(), ite = neighbors26.end(); it != ite; ++it) {
			if (skeletonPoints.find(*it) != skeletonPoints.end())
				cpt++;
		}
		bool add = (cpt <= 2);
		add = true;
		if (centerOfMass != Z3i::Point() && add) {
			markConnectedComponent3D(setVolumeWeighted, connectedComponent3D, i);
			skeletonPoints.insert(centerOfMass);		    
			viewer << CustomColors3D(Color::Red, Color::Red) << centerOfMass;
			viewer << Viewer3D<>::updateDisplay;
			qApp->processEvents();
		}	   

		previousNormal = normal;
		previousCenter = centerOfMass;


		auto pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
				return (wpc->myPoint == (centerOfMass + normal));
			});
		currentPoint = (*pointInWeightedSet);
		if (pointInWeightedSet == setVolumeWeighted.end() || currentPoint->myProcessed) {
			pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
					return (wpc->myPoint == (centerOfMass - normal));
				});
			currentPoint = (*pointInWeightedSet);
			if (pointInWeightedSet == setVolumeWeighted.end() || currentPoint->myProcessed) {
				auto pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
						return (!wpc->myProcessed);
					});
				if (pointInWeightedSet != setVolumeWeighted.end()) {
					currentPoint = (*pointInWeightedSet);
				}
				else break;
			}
		}
		numberLeft = count_if(setVolumeWeighted.begin(), setVolumeWeighted.end(),
							  [&](WeightedPointCount* wpc) {
								  return (!wpc->myProcessed);
							  });
  		i++;
	}
	trace.endBlock();
		
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
	  	viewer << CustomColors3D(Color(0,0,120,10), Color(0,0,120,10)) << (*it)->myPoint;
	}
		
	Image outImage(volume.domain());
	DGtal::imageFromRangeAndValue(skeletonPoints.begin(), skeletonPoints.end(), outImage, 10);
	VolWriter<Image>::exportVol(outFilename, outImage);
	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();
	return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
