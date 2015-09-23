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
#include "geometry/MedialAxis.h"
#include "geometry/ImageUtil.h"
#include "surface/SurfaceUtils.h"
#include "surface/Morphomaths.h"
#include "clustering/diana.hpp"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "geometry/WeightedPointCount.h"
#include "Statistics.h"

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


template <typename WeightedPoint>
set<WeightedPoint*> computePointsBelowPlane(set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume, const Z3i::RealPoint& normal, const Z3i::Point& center, double distanceMax = numeric_limits<double>::max()) {
	set<WeightedPoint*> inferiorPoints;
	double d = -(-normal[0] * center[0] - normal[1] * center[1] - normal[2] * center[2]);
	for (auto it = volume.begin(), ite = volume.end(); it!=ite; ++it) {
		double valueToCheckForPlane = (*it)->myPoint[0] * normal[0] + (*it)->myPoint[1] * normal[1] + (*it)->myPoint[2] * normal[2];
		if (valueToCheckForPlane < d && Z3i::l2Metric((*it)->myPoint, center) <= distanceMax)
			inferiorPoints.insert((*it));
	}
	return inferiorPoints;
}

template <typename WeightedPoint>
void markDifferenceBetweenPlanes(set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume,
											const Z3i::RealPoint& normalPrevious, const Z3i::Point& centerPrevious,
								 const Z3i::RealPoint& normalNext, const Z3i::Point& centerNext,
								 double distanceMax = numeric_limits<double>::max()) {
	Z3i::RealPoint next = normalNext;
	if (normalPrevious.dot(normalNext) < 0) {
		next  = -normalNext;
 	}
	set<WeightedPoint*> inferiorPointsPrevious = computePointsBelowPlane(volume, normalPrevious, centerPrevious, distanceMax);
	set<WeightedPoint*> inferiorPointsNext = computePointsBelowPlane(volume, next, centerNext, distanceMax);
		
	set<WeightedPoint*> difference;
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
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		if (difference.find(*it) != difference.end()) {
			(*it)->myProcessed = true;
		}
	}
}

template <typename Adapter>
Z3i::DigitalSet project2DSetIn3D(const Z2i::DigitalSet& set, const Z3i::Domain& domain3D, const Adapter& adapter) {
	Z3i::DigitalSet set3D(domain3D);
	for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
		set3D.insert(adapter(*it));
	}
	return set3D;
}



bool computePlaneWithGaussianMap(const Z3i::DigitalSet& setSurface, Z3i::DigitalSet& connectedComponent3D, const map<Z3i::Point, Z3i::RealPoint>& normalToFacet, Z3i::RealPoint& normal, const Z3i::Point& currentPoint) {
	typedef Z3i::Object26_6 ObjectType;

	ObjectType object(Z3i::dt26_6, setSurface);
	vector<Z3i::Point> plane;
	unsigned int previousSize = 0;
	int i = 0;

	do {
		previousSize = plane.size();

		auto neighbors = object.neighborhood(currentPoint);
		double d = -(-normal[0] * currentPoint[0] - normal[1] * currentPoint[1] - normal[2] * currentPoint[2]);
		//Naive plane (26 connexity)
		double omega = max(abs(normal[0]), max(abs(normal[1]), abs(normal[2])));
		Z3i::DigitalSet connectedComponent = VCMUtil::extractConnectedComponent3D(setSurface.domain(), setSurface, normal, currentPoint, d, omega);
		for (auto it = connectedComponent.begin(), ite = connectedComponent.end(); it != ite; ++it) {
			if (find(plane.begin(), plane.end(), *it) == plane.end())
				plane.push_back(*it);
		}
		vector<Z3i::RealPoint> fittingPlane;
		for (auto it = plane.begin(), ite = plane.end(); it != ite; ++it) {
			fittingPlane.push_back(normalToFacet.at(*it));
		}
		normal = Statistics::computeNormalFromLinearRegression<Z3i::RealPoint>(fittingPlane);
		i++;
	}
	while (plane.size() != previousSize);

	vector<vector<Z3i::Point> > clusters(plane.size(), vector<Z3i::Point>());
	for (unsigned int i = 0; i < plane.size(); i++) {
		clusters[i].push_back(plane[i]);
	}
	Agnes::mainLoop<Z3i::Point>(clusters);
	double minDistance = numeric_limits<double>::max();
	int index = -1;
	for (unsigned int i = 0; i < clusters.size(); i++) {
		double sum = 0;
		for (auto it = clusters[i].begin(), ite = clusters[i].end(); it != ite; ++it) {
			sum += Z3i::l2Metric(*it, currentPoint);
		}
		sum /= clusters[i].size();
		if (sum < minDistance) {
			minDistance = sum;
			index = i;
		}
	}
	if (index != -1) 
	    plane = clusters[index];
	for (auto it = plane.begin(), ite = plane.end(); it != ite; ++it) {
		connectedComponent3D.insert(*it);
	}
	return (clusters.size() == 1);
}



template <typename Image>
bool isBranchingPart(const Image& image) {
	typedef typename Image::Domain Domain;
	typedef typename Image::Value Scalar;
	
	Domain domain = image.domain();
	Z2i::Point lower = domain.lowerBound(), upper = domain.upperBound();

	int width = upper[0] - lower[0];
	int height = upper[1] - lower[1];

	bool ishape = false, yshape = false;
	for (int i = 0; i < height; i++) {
		int cpt = 0, numberOfBackgroundClusters = 0;
		bool previousZero = false;
		for (int j = 0; j < width; j++) {
			Z2i::Point point(j, i);
			if (domain.isInside(point)) {
				Scalar value = image(point);
			    int valueInt = (int)(value);
				if (valueInt == 0 && previousZero) cpt++;
				else if (valueInt > 0)
				{
					if (cpt >= 0.1*width)
						numberOfBackgroundClusters++;
					cpt = 0;
				}

				previousZero = (valueInt == 0);
			}
		}
		if (numberOfBackgroundClusters == 1)
			ishape = true;
		else if (numberOfBackgroundClusters == 2)
			yshape = true;
		
	}
	return (ishape && yshape);
}

template <typename Image>
bool isBranchingPartConnectedBased(const Image& image) {
	typedef Z2i::Object8_4 ObjectType;
    typedef DGtal::functors::NotPointPredicate<Image> BackgroundPredicate;
	
	Z2i::Domain domain = image.domain();
	Z2i::DigitalSet points(domain);
	BackgroundPredicate bgPredicate(image);
	for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
		if (bgPredicate(*it) > 0)
			points.insert(*it);
	}
	
	ObjectType objectImage(Z2i::dt8_4, points);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);
	return (nbConnectedComponents == 3);
}

template <typename ImageAdapterExtractor, typename Image, typename VCM, typename KernelFunction>
bool isInABranch(const Image& volume, const VCM& vcm, const KernelFunction& chi,
				 const Z3i::Point& point, double radius, int &cptPlane) {
	
	typedef ImageSelector<Z2i::Domain, bool>::Type Image2D;
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type ImageOut;
	
	DGtal::functors::Identity idV;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-radius, -radius, -radius), volume.domain().upperBound() + Z3i::Point(radius, radius, radius));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
									  DGtal::Z2i::Point(radius, radius));
	Z3i::RealPoint longitudinalNormal = computeNormalFromVCM(point, vcm, chi, 1); 
	Z3i::RealPoint sagittalNormal = computeNormalFromVCM(point, vcm, chi, 2);

	DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedderLongitudinal(domain3Dyup, point, longitudinalNormal, radius, domain3Dyup.lowerBound());
	DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedderSagittal(domain3Dyup, point, sagittalNormal, radius, domain3Dyup.lowerBound());
	ImageAdapterExtractor extractedImageLongitudinal(volume, domainImage2D, embedderLongitudinal, idV);
	ImageAdapterExtractor extractedImageSagittal(volume, domainImage2D, embedderSagittal, idV);
	Image2D processImageLongitudinal = ImageUtil::convertImage<Image2D>(extractedImageLongitudinal);
	Image2D processImageSagittal = ImageUtil::convertImage<Image2D>(extractedImageSagittal);

	//Morphological opening
	processImageLongitudinal = Morphomaths::open(processImageLongitudinal, 1);
	processImageSagittal = Morphomaths::open(processImageSagittal, 1);
	
	extractedImageSagittal.setDefaultValue(0);
	extractedImageLongitudinal.setDefaultValue(0);
	bool isBranchingLongitudinal = isBranchingPartConnectedBased(processImageLongitudinal);
	bool isBranchingSagittal = isBranchingPartConnectedBased(processImageSagittal);

	if (isBranchingSagittal) {
		string outName = "/home/florent/trash/slice_" + std::to_string(cptPlane) + ".pgm";
		ImageOut outSagittal = ImageUtil::convertImage<ImageOut>(processImageSagittal, 255);
		PGMWriter<ImageOut>::exportPGM(outName, outSagittal);
		cptPlane++;
	} if (isBranchingLongitudinal) {
		string outName = "/home/florent/trash/slice_" + std::to_string(cptPlane) + ".pgm";
		ImageOut outLongitudinal = ImageUtil::convertImage<ImageOut>(processImageLongitudinal, 255);
		PGMWriter<ImageOut>::exportPGM(outName, outLongitudinal);
		cptPlane++;

	}
	
	return (isBranchingSagittal || isBranchingLongitudinal);
	
}

template <typename VCM, typename KernelFunction, typename Container, typename DT>
void connectDisconnectedComponents(Z3i::DigitalSet& skeletonPoints, const DT& dt, double delta,
								   VCM& vcm,  KernelFunction& chi,
								   const Z3i::DigitalSet& setVolume, Container& setVolumeWeighted) {
	typedef Z3i::Object26_6 ObjectType;
	typedef WeightedPointCount<Z3i::Point> WeightedPointCount;

	set<Z3i::Point> pointsToProcess;
	Z3i::Domain domainVolume = setVolume.domain();
	ObjectType objectImage(Z3i::dt26_6, skeletonPoints);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);
	trace.info() << nbConnectedComponents << endl;

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
		Z3i::RealPoint directionVector = (belongingToReference - belongingToCurrentObject).getNormalized();
		if (directionVector.dot(normal) < 0)
			normal = -normal;
		map<Z3i::Point, int> accumulator;
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
				accumulator[referencePoint]++;
			}
			scalar++;	
			previousProjection = projection;
		}

		std::vector<pair<Z3i::Point, int>> pairs;
		for (auto itr = accumulator.begin(); itr != accumulator.end(); ++itr)
			pairs.push_back(*itr);

		sort(pairs.begin(), pairs.end(), [&](const pair<Z3i::Point, int>& a, const pair<Z3i::Point, int>& b)
			 {
				 return a.second > b.second;
			 }
			);
		belongingToReference = (pairs.begin())->first;
		trace.info() << (pairs.begin())->second;
		for (auto it = minimizingObjectToReference->pointSet().begin(),
				 ite = minimizingObjectToReference->pointSet().end();
			 it != ite; ++it) {
			reference.pointSet().insert(*it);
		}
		
		vector<Z3i::Point> points = PointUtil::linkTwoPoints(belongingToCurrentObject, belongingToReference);
		for (auto it = points.begin(), ite = points.end(); it!=ite; ++it) {
			pointsToProcess.insert(*it);
			reference.pointSet().insert(*it);
		}
		objects.erase(minimizingObjectToReference);
	}

	
	
	double radius = 0;
	Z3i::RealPoint normal, previousNormal;
	int i = 0;
	for (auto it = pointsToProcess.begin(), ite = pointsToProcess.end();
		 it != ite; ++it) {
		trace.progressBar(i, pointsToProcess.size());
		radius = dt(*it) + delta;
		vcm.updateProximityStructure(radius, setVolume.begin(), setVolume.end());
		chi = KernelFunction( 1.0, radius);
		Z3i::DigitalSet connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted,
																	*it, normal, 0, previousNormal, radius);
		Z3i::Point centerOfMass = Statistics::extractCenterOfMass3D(connectedComponent3D);
				
		if (centerOfMass != Z3i::Point()) {
			reference.pointSet().insert(centerOfMass);
			skeletonPoints.insert(centerOfMass);
		}
		i++;
	}
	trace.endBlock();
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
		("input,i", po::value<std::string>(), "vol file (corresponding volume)")
		("output,o", po::value<std::string>(), "output skeleton filename")
		("skeleton,s", po::bool_switch()->default_value(false), "vol file (medial axis)")
		("delta,d", po::value<int>()->default_value(1), "delta for ball radius")
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
	int delta = vm["delta"].as<int>();
	int cptPlane = 0;
	
	Image volume = VolReader<Image>::importVol(inputFilename);
	Z3i::Domain domainVolume = volume.domain();
	Z3i::DigitalSet setVolume(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume, 
												  thresholdMin-1, thresholdMax);
	set<WeightedPointCount*, WeightedPointCountComparator<WeightedPointCount>> setVolumeWeighted;
	Image volumeBinary(volume.domain());
	for (auto it = volume.domain().begin(), ite = volume.domain().end(); it != ite; ++it) {
		if (volume(*it) >= thresholdMin && volume(*it) <= thresholdMax) 
			volumeBinary.setValue(*it, 255);
		else
			volumeBinary.setValue(*it, 0);
	}

	vector<Point> vPoints;
	Z3i::DigitalSet skeletonPoints(domainVolume);
	typedef DGtal::functors::IntervalForegroundPredicate<Image> Binarizer; 
	Binarizer binarizer(volume, thresholdMin-1, thresholdMax);
	typedef DGtal::DistanceTransformation<Space, Binarizer, Z3i::L2Metric> DTL2;
	
	DTL2 dt(&volume.domain(), &binarizer, &Z3i::l2Metric);
		

	DTL2::Domain domainDT = dt.domain();
	for (auto it = domainDT.begin(), ite = domainDT.end(); it != ite; ++it) {
		double value = dt(*it);
		if (value > 0) {
			setVolumeWeighted.insert(new WeightedPointCount(*it, value));
			checkPointForMedialAxis(dt, vPoints, *it);
		}		
	}
	trace.info() << vPoints.size() << endl;
	WeightedPointCount* currentPoint = *setVolumeWeighted.begin();
  
	Metric l2;
	VCM vcm( R, ceil( r ), l2, true );
	vcm.init( setVolume.begin(), setVolume.end() );
	Domain domain = vcm.domain();
	KernelFunction chi( 1.0, r );

	
	Matrix vcm_r, evec;
	RealVector eval;

	Viewer3D<> viewer;
	viewer.show();
   
 
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );

	int i = 0;
	int numberLeft = setVolumeWeighted.size();
	
	Z3i::RealPoint normal;
	Z3i::RealPoint previousNormal;

	Z3i::DigitalSet connectedComponent3D(domainVolume);
	Z3i::RealPoint realCenter;
	Z3i::Point centerOfMass;
	Z3i::Point previousCenter;
	
	trace.beginBlock("Computing skeleton");
	//Main loop to compute skeleton (stop when no vol points left to process)
	while (numberLeft > 0)
	{
		trace.progressBar((setVolume.size() - numberLeft), setVolume.size());		
		currentPoint->myProcessed = true;
		double radius = r;

		//Distance transform value for VCM radius
		if (vm.count("skeleton")) {
			Point closestPointToCurrent = *min_element(vPoints.begin(), vPoints.end(), [&](const Point& one, const Point& two) {
					return Z3i::l2Metric(one, currentPoint->myPoint) < Z3i::l2Metric(two, currentPoint->myPoint);
				});
			radius = dt(closestPointToCurrent);
			radius += delta;
			if (radius > 0) {
				vcm.updateProximityStructure(radius, setVolume.begin(), setVolume.end());
				chi = KernelFunction( 1.0, radius);
			}
		}
	
		// Compute discrete plane
		connectedComponent3D = VCMUtil::computeDiscretePlane(vcm, chi, domainVolume, setVolumeWeighted, currentPoint->myPoint, normal,
																	0, previousNormal, radius);
	    realCenter = Statistics::extractCenterOfMass3D(connectedComponent3D);
	

		//Branching detection
		//bool inABranch = isInABranch<ImageAdapterExtractor>(volumeBinary, vcm, chi, currentPoint->myPoint, (radius)*5, cptPlane);
		// if (inABranch) {
		// 	viewer << CustomColors3D(Color::Blue, Color::Blue) << currentPoint->myPoint;
		// }
		
		//Center of mass computation
		if (realCenter != Z3i::RealPoint()) {
		    bool add = VCMUtil::markConnectedComponent3D(setVolumeWeighted, connectedComponent3D, i);
			if (add) {
				centerOfMass = extractNearestNeighborInSetFromPoint(connectedComponent3D, realCenter);
				previousNormal = normal;
				previousCenter = centerOfMass;
				skeletonPoints.insert(centerOfMass);
				viewer << CustomColors3D(Color::Red, Color::Red) << centerOfMass;
				viewer << Viewer3D<>::updateDisplay;
				qApp->processEvents();
			}
		}
		
		//Go to next point according to normal OR to max value in DT
		auto pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
				return (wpc->myPoint == centerOfMass);
			});
		int scalar = 1;
		Point current = centerOfMass;
	    auto newPoint = setVolumeWeighted.begin();
		if (pointInWeightedSet != setVolumeWeighted.end()) {
			while (current == centerOfMass ||
				   connectedComponent3D.find(current) != connectedComponent3D.end()) {
				current = centerOfMass + normal * scalar;
				scalar++;
			}
			newPoint = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
					return (wpc->myPoint == current);
				});
			if (newPoint == setVolumeWeighted.end() || (*newPoint)->myProcessed) {
				scalar = 1;
				current = centerOfMass;
				if (pointInWeightedSet != setVolumeWeighted.end()) {
					while (current == centerOfMass ||
						connectedComponent3D.find(current) != connectedComponent3D.end()) {
						current = centerOfMass - normal * scalar;
						scalar++;
					}
					newPoint = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
							return (wpc->myPoint == current);
						});
					if (newPoint == setVolumeWeighted.end() || (*newPoint)->myProcessed) {
						previousNormal = Z3i::RealPoint();
						pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
								return (!wpc->myProcessed);
							});
						if (pointInWeightedSet != setVolumeWeighted.end()) {
							newPoint = pointInWeightedSet;
						}
					}
				}
			}
		} else {
			previousNormal = Z3i::RealPoint();
			pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
					return (!wpc->myProcessed);
				});
			if (pointInWeightedSet != setVolumeWeighted.end()) {
				newPoint = pointInWeightedSet;
			}
		}		

		currentPoint = (*newPoint);
		
		numberLeft = count_if(setVolumeWeighted.begin(), setVolumeWeighted.end(),
							  [&](WeightedPointCount* wpc) {
								  return (!wpc->myProcessed);
							  });
  		i++;
	}
	trace.endBlock();

	viewer << Viewer3D<>::updateDisplay;
	qApp->processEvents();

	//Second pass
	for (auto it = setVolumeWeighted.begin(), ite = setVolumeWeighted.end(); it != ite; ++it) {
		(*it)->myProcessed = false;
	}
	connectDisconnectedComponents(skeletonPoints, dt, delta, vcm, chi, setVolume, setVolumeWeighted);
   
	//Displaying
	for (auto it = skeletonPoints.begin(), ite = skeletonPoints.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color::Green, Color::Green) << *it;
	}

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
