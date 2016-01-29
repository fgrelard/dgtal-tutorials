#ifndef VCM_UTIL_H
#define VCM_UTIL_H

#include "DGtal/base/Common.h"

#include "DGtal/helpers/StdDefs.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "WeightedPointCountComparator.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/WeightedPointCount.h"
#include "geometry/PointUtil.h"

#include <set>
#include <vector>

namespace VCMUtil {
	template <typename VCM, typename KernelFunction>
	DGtal::Z3i::RealPoint computeNormalFromVCM(const DGtal::Z3i::Point& currentPoint, const VCM& vcm, const KernelFunction& chi, int coordinate, const DGtal::Z3i::RealVector& dirVector = DGtal::Z3i::RealVector());

	template <typename VCM, typename KernelFunction>
	DGtal::Z3i::RealPoint computeEigenValuesFromVCM(const DGtal::Z3i::Point& currentPoint, const VCM& vcm, const KernelFunction& chi);


	bool planeContains(const DGtal::Z3i::Point& point, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& current);
	bool abovePlane(const DGtal::Z3i::Point& point, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& current);

	
	template <typename Domain, typename WeightedPoint>
	void extractConnectedComponent3D(DGtal::Z3i::DigitalSet& intersection, const Domain & domain, const std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint> >& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double& distanceMax);

	template <typename Image>
	DGtal::Z2i::DigitalSet extractConnectedComponent(const Image& image, const DGtal::Z2i::Point& referencePoint, int thresholdMin,
														  int thresholdMax);

	template <typename Domain>
	DGtal::Z3i::DigitalSet extractConnectedComponent3D(const Domain & domain, const DGtal::Z3i::DigitalSet& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double distanceMax = std::numeric_limits<double>::max());

	bool isRadiusMaximal(const DGtal::Z3i::DigitalSet& intersection, const DGtal::Z3i::Point& referencePoint, double currentDistance, const DGtal::Z3i::RealVector& dirVector);
	
	template <typename WeightedPoint>
	bool markConnectedComponent3D(std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint>>& volume, const DGtal::Z3i::DigitalSet& intersection, int label);

	template <typename VCM, typename KernelFunction, typename Domain, typename Container>
	DGtal::Z3i::DigitalSet computeDiscretePlane(VCM& vcm, KernelFunction& chi, const Domain& domainVolume, const Container& setVolumeWeighted, const DGtal::Z3i::Point& point, DGtal::Z3i::RealPoint& normal, int coordinate, double& radius, double distanceMax, bool dilate=true, const DGtal::Z3i::RealVector dirVector=DGtal::Z3i::RealVector());

	template <typename WeightedPointCount, typename Container>
	void trackNextPoint(WeightedPointCount* &currentPoint, const Container& setVolumeWeighted,
						const DGtal::Z3i::DigitalSet& connectedComponent3D,
						const DGtal::Z3i::Point& centerOfMass, const DGtal::Z3i::RealPoint& normal);

	double computeCurvatureJunction(const DGtal::Z3i::RealPoint& lambda);

	template <typename VCM>
	double radiusAtJunction(const VCM& vcm, const DGtal::Z3i::RealPoint& branchPoint, double radius);
}

template <typename VCM, typename KernelFunction>
DGtal::Z3i::RealPoint VCMUtil::computeNormalFromVCM(const DGtal::Z3i::Point& currentPoint, const VCM& vcm, const KernelFunction& chi, int coordinate, const DGtal::Z3i::RealVector& dirVector) {
	
	typedef DGtal::EigenDecomposition<3,double> LinearAlgebraTool;
	LinearAlgebraTool::Matrix vcm_r, evec;
	DGtal::Z3i::RealVector eval;
// Compute VCM and diagonalize it.
	if (dirVector == DGtal::Z3i::RealVector())
		vcm_r = vcm.measure( chi, currentPoint);
	else
		vcm_r = vcm.measureJunction( dirVector, chi, currentPoint);
	LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
	// // Display normal
	DGtal::Z3i::RealVector normal = evec.column(coordinate);
	return normal;
}

template <typename VCM, typename KernelFunction>
DGtal::Z3i::RealPoint VCMUtil::computeEigenValuesFromVCM(const DGtal::Z3i::Point& currentPoint, const VCM& vcm, const KernelFunction& chi) {
	typedef DGtal::EigenDecomposition<3,double> LinearAlgebraTool;
	LinearAlgebraTool::Matrix vcm_r, evec;
	DGtal::Z3i::RealVector eval;
// Compute VCM and diagonalize it.
	vcm_r = vcm.measure( chi, currentPoint );
	LinearAlgebraTool::getEigenDecomposition( vcm_r, evec, eval );
    
	return eval;
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

bool VCMUtil::planeContains(const DGtal::Z3i::Point& point, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& current) {
	double d = -(-normal[0] * current[0] - normal[1] * current[1] - normal[2] * current[2]);
	//Naive plane (26 connexity)
	double omega = std::max(abs(normal[0]), std::max(abs(normal[1]), abs(normal[2])));
	double valueToCheckForPlane = point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2];
	if (valueToCheckForPlane >= d && valueToCheckForPlane < d + omega)
		return true;
	return false;
}

bool VCMUtil::abovePlane(const DGtal::Z3i::Point& point, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& current) {
	double d = -(-normal[0] * current[0] - normal[1] * current[1] - normal[2] * current[2]);
	//Naive plane (26 connexity)
	double valueToCheckForPlane = point[0] * normal[0] + point[1] * normal[1] + point[2] * normal[2];
	if (valueToCheckForPlane >= d)
		return true;
	return false;
}

template <typename Domain, typename WeightedPoint>
void VCMUtil::extractConnectedComponent3D(DGtal::Z3i::DigitalSet& intersection, const Domain & domain, const std::set<WeightedPoint*, WeightedPointCountComparator<WeightedPoint> >& volume, const DGtal::Z3i::RealPoint& normal, const DGtal::Z3i::Point& referencePoint, double d, double omega, double& distanceMax) {
	typedef DGtal::Z3i::Object26_6 ObjectType;

	DGtal::Z3i::DigitalSet aSet(domain);
	for (auto it = volume.begin(), ite = volume.end(); it != ite; ++it) {
		double valueToCheckForPlane = (*it)->myPoint[0] * normal[0] + (*it)->myPoint[1] * normal[1] + (*it)->myPoint[2] * normal[2];
		if (valueToCheckForPlane >= d && valueToCheckForPlane < d + omega) {
			aSet.insert((*it)->myPoint);
		}
	}

	ObjectType objectIntersection(DGtal::Z3i::dt26_6, aSet);
	std::vector<ObjectType> objects;
	std::back_insert_iterator<std::vector<ObjectType>> inserter(objects);
	unsigned int nbConnectedComponents = objectIntersection.writeComponents(inserter);
	if (nbConnectedComponents > 1) {
		for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
			DGtal::Z3i::DigitalSet ccSet = it->pointSet();		    
			if (ccSet.find(referencePoint) != ccSet.end()) {
				intersection = ccSet;
			}
		}
	} else {
		intersection = aSet;
	}

}

bool VCMUtil::isRadiusMaximal(const DGtal::Z3i::DigitalSet& intersection, const DGtal::Z3i::Point& referencePoint, double currentDistance, const DGtal::Z3i::RealVector& dirVector) {
	bool alright = true;
	if (dirVector == DGtal::Z3i::RealVector()) {
		for (auto it = intersection.begin(), ite = intersection.end(); it != ite; ++it) {
			if (DGtal::Z3i::l2Metric(*it, referencePoint) > currentDistance) {
				alright = false;
			}
		}
	} else {
		DGtal::Z3i::Point current = referencePoint;
		int scalar = 1;
		while (intersection.find(current) != intersection.end()) {
			current = referencePoint + dirVector * scalar;
			scalar++;
		}
		double distance = DGtal::Z3i::l2Metric(current, referencePoint);
		alright = (distance <= currentDistance);
	}
   
	return alright;
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
DGtal::Z3i::DigitalSet VCMUtil::computeDiscretePlane(VCM& vcm, KernelFunction& chi,
													 const Domain& domainVolume, const Container& setVolumeWeighted,
													 const DGtal::Z3i::Point& point, DGtal::Z3i::RealPoint& normal, int coordinate, double& radius,
													 double distanceMax, bool dilate, const DGtal::Z3i::RealVector dirVector) {

	bool alright = false;
	DGtal::Z3i::DigitalSet connectedComponent3D(domainVolume);
	double currentDistance = radius;
	do  {
	    currentDistance++;
		chi = KernelFunction(1.0, currentDistance);
		vcm.setMySmallR(currentDistance);
		if (dirVector == DGtal::Z3i::RealVector())
			normal = computeNormalFromVCM(point, vcm, chi, coordinate);
		else
			normal = computeNormalFromVCM(point, vcm, chi, coordinate, dirVector);
		
		double d = -(-normal[0] * point[0] - normal[1] * point[1] - normal[2] * point[2]);
		//Naive plane (26 connexity)
		double omega = std::max(std::abs(normal[0]), std::max(std::abs(normal[1]), std::abs(normal[2])));
		connectedComponent3D = DGtal::Z3i::DigitalSet(domainVolume);
	    extractConnectedComponent3D(connectedComponent3D, domainVolume, setVolumeWeighted, normal, point, d, omega, currentDistance);
		alright = isRadiusMaximal(connectedComponent3D, point, currentDistance, dirVector);		
	} while (!alright && dilate && currentDistance < distanceMax);
	
	radius = currentDistance;

	DGtal::Z3i::DigitalSet discretePlane(domainVolume);
	for (auto it = connectedComponent3D.begin(), ite = connectedComponent3D.end(); it != ite; ++it) {
		if (DGtal::Z3i::l2Metric(*it, point) <= radius)
			discretePlane.insert(*it);
	}
	return discretePlane;
}

template <typename WeightedPointCount, typename Container>
void VCMUtil::trackNextPoint(WeightedPointCount* &currentPoint, const Container& setVolumeWeighted,
					const DGtal::Z3i::DigitalSet& connectedComponent3D,
					const DGtal::Z3i::Point& centerOfMass, const DGtal::Z3i::RealPoint& normal) {
	auto pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
			return (wpc->myPoint == centerOfMass);
		});
	int scalar = 1;
	DGtal::Z3i::Point current = centerOfMass;
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
		pointInWeightedSet = find_if(setVolumeWeighted.begin(), setVolumeWeighted.end(), [&](WeightedPointCount* wpc) {
				return (!wpc->myProcessed);
			});
		if (pointInWeightedSet != setVolumeWeighted.end()) {
			newPoint = pointInWeightedSet;
		}
	}
	currentPoint = (*newPoint);
}

double VCMUtil::computeCurvatureJunction(const DGtal::Z3i::RealPoint& lambda) {
	double ratio = (lambda[0]) / (lambda[0] + lambda[1] + lambda[2]);
	return ratio;
}

template <typename VCM>
double VCMUtil::radiusAtJunction(const VCM& vcm, const DGtal::Z3i::RealPoint& branchPoint, double radius) {
	auto lambda = ((vcm.mapPoint2ChiVCM()).at(branchPoint)).values;
	double ratio = computeCurvatureJunction(lambda);
	double r = radius * (1/(1-ratio*3));
	return r;
}

#endif
