#ifndef INTEGRAL_ADAPTED_H
#define INTEGRAL_ADAPTED_H

#include "DGtal/geometry/surfaces/estimation/IntegralInvariantVolumeEstimator.h"

using namespace DGtal;
using namespace std;

template <typename TKSpace, typename TPointPredicate, typename TCovarianceMatrixFunctor>
class IntegralInvariantVolumeEstimatorAdapted : public IntegralInvariantVolumeEstimator<TKSpace, TPointPredicate, TCovarianceMatrixFunctor> {

	using IntegralInvariantVolumeEstimator<TKSpace, TPointPredicate, TCovarianceMatrixFunctor>::IntegralInvariantVolumeEstimator;
	
public:
	template <typename Vertex, typename WeightedSurfel, typename OutputIterator, typename SurfelConstIterator>
	void eval(OutputIterator&, const vector<WeightedSurfel> &, float, float);

	// --- Private methods ---
private:
	template <typename Vertex, typename WeightedSurfel>
	vector<Vertex> findLocation( const vector<WeightedSurfel> &, float, float);

	template <typename WeightedSurfel>
	vector<float> differentValues(const vector<WeightedSurfel> &, float);
};

#include "IntegralInvariantVolumeEstimatorAdapted.ih"

#endif
