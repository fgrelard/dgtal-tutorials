#ifndef SURFACE_UTILS_H
#define SURFACE_UTILS_H

#include <vector>
#include <set>
#include <map>


namespace SurfaceUtils {
	template <typename SurfacePoint, typename KSpace, typename Surfel, typename Point, typename DigitalSet>
	std::vector<SurfacePoint> computeSurfelWeight(const KSpace & ks, const DigitalSet & set3d, const std::set<Surfel>& surfelSet, std::set<Point>& surfacePointSet, std::map<Surfel, Point> & surfaceVoxelSet);
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

#endif
