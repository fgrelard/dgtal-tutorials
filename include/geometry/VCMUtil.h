#ifndef VCM_UTIL_H
#define VCM_UTIL_H

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "WeightedPointCountComparator.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include <set>
#include <vector>

namespace VCMUtil {
	template <typename VCM, typename KernelFunction>
	DGtal::Z3i::RealPoint computeNormalFromVCM(const DGtal::Z3i::Point& currentPoint, VCM& vcm, KernelFunction& chi, int coordinate);

	template <typename Domain, typename WeightedPoint>
	DGtal::Z3i::DigitalSet extractConnectedComponent3D(const Domain & domain, const std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint> >& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double distanceMax);

	template <typename Image>
	DGtal::Z2i::DigitalSet extractConnectedComponent(const Image& image, const DGtal::Z2i::Point& referencePoint, int thresholdMin,
														  int thresholdMax);

	template <typename Domain>
	DGtal::Z3i::DigitalSet extractConnectedComponent3D(const Domain & domain, const DGtal::Z3i::DigitalSet& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double distanceMax = std::numeric_limits<double>::max());
	
	template <typename WeightedPoint>
	bool markConnectedComponent3D(std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume, const DGtal::Z3i::DigitalSet& intersection, int label);

	template <typename VCM, typename KernelFunction, typename Domain, typename Container>
	DGtal::Z3i::DigitalSet computeDiscretePlane(const VCM& vcm, const KernelFunction& chi, const Domain& domainVolume, const Container& setVolumeWeighted, const DGtal::Z3i::Point& point, DGtal::Z3i::RealPoint& normal, int coordinate, const DGtal::Z3i::RealPoint previousNormal, double& radius);
	
}

template <typename VCM, typename KernelFunction>
DGtal::Z3i::RealPoint VCMUtil::computeNormalFromVCM(const DGtal::Z3i::Point& currentPoint, VCM& vcm, KernelFunction& chi, int coordinate) {
	
	typedef DGtal::EigenDecomposition<3,double> LinearAlgebraTool;
	LinearAlgebraTool::Matrix vcm_r, evec;
	DGtal::Z3i::RealVector eval;
// Compute VCM and diagonalize it.
	vcm_r = vcm.measure( chi, currentPoint );
	LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
	
	// // Display normal
	DGtal::Z3i::RealVector normal = evec.column(coordinate);
	return normal;
}

template <typename Domain>
DGtal::Z3i::DigitalSet VCMUtil::extractConnectedComponent3D(const Domain & domain, const DGtal::Z3i::DigitalSet& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double distanceMax = std::numeric_limits<double>::max()) {
	typedef DGtal::Z3i::Object26_6 ObjectType;

	DGtal::Z3i::DigitalSet intersection(domain);
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		double valueToCheckForPlane = (*it)[0] * normal[0] + (*it)[1] * normal[1] + (*it)[2] * normal[2];
		if (valueToCheckForPlane >= d && valueToCheckForPlane < d + omega &&
			DGtal::Z3i::l2Metric((*it), referencePoint) <= distanceMax) {
			intersection.insert((*it));
		}
	}

	ObjectType objectIntersection(DGtal::Z3i::dt26_6, intersection);
	std::vector<ObjectType> objects;
	std::back_insert_iterator<std::vector<ObjectType>> inserter(objects);
	unsigned int nbConnectedComponents = objectIntersection.writeComponents(inserter);
	DGtal::Z3i::DigitalSet connectedComponent = intersection;
	double min = std::numeric_limits<double>::max();
	if (nbConnectedComponents > 1) {
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			double sum = 0;
			DGtal::Z3i::DigitalSet ccSet = it->pointSet();
		    for (auto it = ccSet.begin(), ite = ccSet.end(); it != ite; ++it) {
				sum += DGtal::Z3i::l2Metric(*it, referencePoint);
			}
			sum /= ccSet.size();
			if (sum < min) {
				min = sum;
				connectedComponent = ccSet;
			}
		}
	}
	return connectedComponent;
}

template <typename Domain, typename WeightedPoint>
DGtal::Z3i::DigitalSet VCMUtil::extractConnectedComponent3D(const Domain & domain, const std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint> >& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double distanceMax) {
	typedef DGtal::Z3i::Object26_6 ObjectType;

	DGtal::Z3i::DigitalSet intersection(domain);
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		double valueToCheckForPlane = (*it)->myPoint[0] * normal[0] + (*it)->myPoint[1] * normal[1] + (*it)->myPoint[2] * normal[2];
		if (valueToCheckForPlane >= d && valueToCheckForPlane < d + omega && DGtal::Z3i::l2Metric((*it)->myPoint, referencePoint) <= distanceMax) {
			intersection.insert((*it)->myPoint);
		}
	}

	ObjectType objectIntersection(DGtal::Z3i::dt26_6, intersection);
	std::vector<ObjectType> objects;
	std::back_insert_iterator<std::vector<ObjectType>> inserter(objects);
	unsigned int nbConnectedComponents = objectIntersection.writeComponents(inserter);
	if (nbConnectedComponents > 1) {
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			DGtal::Z3i::DigitalSet ccSet = it->pointSet();
			if (ccSet.find(referencePoint) != ccSet.end()) {
				
				return ccSet;
			}
		}
		return DGtal::Z3i::DigitalSet(domain);
	}
   
	return intersection;
}



template <typename Image>
DGtal::Z2i::DigitalSet VCMUtil::extractConnectedComponent(const Image& image, const DGtal::Z2i::Point& referencePoint, int thresholdMin,
										  int thresholdMax) {
	typedef DGtal::ImageSelector<DGtal::Z2i::Domain, unsigned char>::Type Image2D;
	typedef DGtal::Z2i::Object8_4 ObjectType;
	
	//Extracting connected components
	DGtal::Z2i::DigitalSet points2D(image.domain());
	DGtal::SetFromImage<DGtal::Z2i::DigitalSet>::append<Image> (points2D, image, 
												  thresholdMin-1, thresholdMax);
	ObjectType objectImage2D(DGtal::Z2i::dt8_4, points2D);
	std::vector<ObjectType> objects;
	std::back_insert_iterator< std::vector< ObjectType > > inserter( objects );
	unsigned int nbConnectedComponents = objectImage2D.writeComponents(inserter);
	Image2D image2D(image.domain());
	if (nbConnectedComponents > 1) {
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			DGtal::Z2i::DigitalSet ccSet = it->pointSet();
			if (ccSet.find(referencePoint) != ccSet.end()) {
			    return ccSet;
			}
		}
		return DGtal::Z2i::DigitalSet(image.domain());
	}
	return points2D;
}

template <typename WeightedPoint>
bool VCMUtil::markConnectedComponent3D(std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume, const DGtal::Z3i::DigitalSet& intersection, int label) {	

	std::set<int> processedLabels;
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		if (!(*it)->myProcessed && intersection.find((*it)->myPoint) != intersection.end()) {
			(*it)->myProcessed = true;
			if ((*it)->myCount != 1)
				(*it)->myCount = label;
		}
		else if ((*it)->myProcessed  && intersection.find((*it)->myPoint) != intersection.end()) {
			processedLabels.insert((*it)->myCount);
		}
	}
	return (processedLabels.size() <= 2);	
}


template <typename VCM, typename KernelFunction, typename Domain, typename Container>
DGtal::Z3i::DigitalSet VCMUtil::computeDiscretePlane(const VCM& vcm, const KernelFunction& chi,
							   const Domain& domainVolume, const Container& setVolumeWeighted,
									 const DGtal::Z3i::Point& point, DGtal::Z3i::RealPoint& normal, int coordinate,
									 const DGtal::Z3i::RealPoint previousNormal, double& radius) {

    normal = computeNormalFromVCM(point, vcm, chi, coordinate);
	if (normal.dot(previousNormal) < 0)
		normal = -normal;
		
		
	double d = -(-normal[0] * point[0] - normal[1] * point[1] - normal[2] * point[2]);
	//Naive plane (26 connexity)
	double omega = std::max(abs(normal[0]), std::max(abs(normal[1]), abs(normal[2])));
	DGtal::Z3i::DigitalSet connectedComponent3D = extractConnectedComponent3D(domainVolume, setVolumeWeighted, normal, point, d, omega, radius);

   
	return connectedComponent3D;
}

#endif
