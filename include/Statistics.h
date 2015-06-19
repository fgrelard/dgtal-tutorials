#ifndef __STATISTICS__
#define __STATISTICS
#include <vector>
#include <algorithm>

namespace Statistics {
	template <typename T>
	double mean(const std::vector<T>& aVector) {
		return std::accumulate( aVector.begin(), aVector.end(), 0.0f )/ aVector.size();
	}

	template <typename T>
	double stddev(const std::vector<T>& aVector) {
		std::vector<double> zero_mean( aVector );
		transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( std::minus<double>(), mean(aVector) ) );

		double deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
		deviation = sqrt( deviation / ( aVector.size() - 1 ) );
		return deviation;
	}
}
#endif
