#ifndef DSS_UTIL_H
#define DSS_UTIL_H

#include <vector>
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/base/Common.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/helpers/StdDefs.h"

namespace DSSUtil {

	template <typename Segmentation>
	std::vector<typename Segmentation::SegmentComputer> computeDSSPassingThrough(const DGtal::Z3i::Point& point, const Segmentation& segmentation);

    template <typename DSS>
	std::vector<double> extractLengths(const std::vector<DSS>& vDSS);
	 
}

template <typename Segmentation>
std::vector<typename Segmentation::SegmentComputer> DSSUtil::computeDSSPassingThrough(const DGtal::Z3i::Point& point, const Segmentation& segmentation) {
	typedef typename Segmentation::SegmentComputer DSS;
	std::vector<DSS> vDSS;
	for (auto it = segmentation.begin(), ite = segmentation.end(); it != ite; ++it) {
		for (auto itP = it->begin(), itPe = it->end(); itP != itPe;++itP) {
			if (*itP == point) {
				DSS currentSegment(*it);
				vDSS.push_back(currentSegment);
			}
		}
	}
	return vDSS;

}

template <typename DSS>
std::vector<double> DSSUtil::extractLengths(const std::vector<DSS>& vDSS) {
	std::vector<double> lengths;
	for (auto it = vDSS.begin(), ite = vDSS.end(); it != ite; ++it) {
		DSS currentSegment = *it;
		DGtal::Z3i::RealPoint v = *(--currentSegment.end()) - *currentSegment.begin();
		lengths.push_back(v.norm());
	}
	return lengths;
}

#endif
