#ifndef SURFACE_UTILS_H
#define SURFACE_UTILS_H

#include <vector>
#include <set>
#include <map>
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"

namespace SurfaceUtils {
	template <typename SurfacePoint, typename KSpace, typename Surfel, typename Point, typename DigitalSet>
	std::vector<SurfacePoint> computeSurfelWeight(const KSpace & ks, const DigitalSet & set3d, const std::set<Surfel>& surfelSet, std::set<Point>& surfacePointSet, std::map<Surfel, Point> & surfaceVoxelSet);

	template <typename Image>
	DGtal::Z3i::DigitalSet extractSurfaceVoxels(const Image& volume, int thresholdMin, int thresholdMax);

	template <typename KSpace>
	std::vector<DGtal::Z3i::RealPoint> normalsToSurfel(const typename KSpace::SCell&  surfel);
}

template <typename SurfacePoint, typename KSpace, typename Surfel, typename Point, typename DigitalSet>
std::vector<SurfacePoint> SurfaceUtils::computeSurfelWeight(const KSpace & ks, const DigitalSet & set3d, const std::set<Surfel>& surfelSet, std::set<Point>& surfacePointSet, std::map<Surfel, Point> & surfaceVoxelSet) {
	std::vector<SurfacePoint> weightedSurfaceVector;
	for (auto it = set3d.begin(), itE = set3d.end(); it != itE; ++it) {
		std::vector<Surfel> aSurfelV;
		Surfel current = ks.sSpel(*it);
		int number = 0;
		for (int i = 0; i < 3; i++) {
			auto itSurfel = surfelSet.find(ks.sIncident(current, i, true));
			if (itSurfel != surfelSet.end()) {
				number++;
				aSurfelV.push_back(*itSurfel);
				surfaceVoxelSet[*itSurfel] = *it;
				surfacePointSet.insert(*it);
			}
		}

		for (int i = 0; i < 3; i++) {
			auto itSurfel = surfelSet.find(ks.sIncident(current, i, false));
			if (itSurfel != surfelSet.end()) {
				number++;
				aSurfelV.push_back(*itSurfel);
				surfaceVoxelSet[*itSurfel] = *it;
				surfacePointSet.insert(*it);
			}
		}
		weightedSurfaceVector.push_back({*it, number, aSurfelV});
	}
	return weightedSurfaceVector;
}

template <typename Image>
DGtal::Z3i::DigitalSet SurfaceUtils::extractSurfaceVoxels(const Image& volume, int thresholdMin, int thresholdMax) {

	DGtal::Z3i::DigitalSet surfacePoints(volume.domain());
		
	typedef DGtal::functors::IntervalForegroundPredicate<Image> Binarizer; 
	Binarizer binarizer(volume, thresholdMin-1, thresholdMax);
	typedef DGtal::DistanceTransformation<DGtal::Z3i::Space, Binarizer, DGtal::Z3i::L2Metric> DTL2;
	
	DTL2 dt(&volume.domain(), &binarizer, &DGtal::Z3i::l2Metric);

	for (auto it = volume.domain().begin(), ite = volume.domain().end(); it != ite; ++it) {
		if (dt(*it) == 1)
			surfacePoints.insert(*it);
	}
	return surfacePoints;
}

template <typename KSpace>
std::vector<DGtal::Z3i::RealPoint> normalsToSurfel(const KSpace ks, const typename KSpace::SCell&  surfel) {
	std::vector<DGtal::Z3i::RealPoint> normals;
	for (auto it = ks.sDirs(surfel); it != 0; ++it) {
		DGtal::Z3i::RealPoint p(0, 0, 0);
		p[*it] = 1;
		normals.push_back(p);
	}
	return normals;
}


#endif
