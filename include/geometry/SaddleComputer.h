#ifndef SADDLE_COMPUTER_H
#define SADDLE_COMPUTER_H

#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "geometry/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "geometry/VCMUtil.h"

namespace SaddleComputer {
	template <typename VCM, typename Domain>
	DGtal::Z3i::DigitalSet computeBranchingPartsWithVCMFeature(const VCM& vcm,
														const Domain& domain, double threshold);

	template <typename VCM>
	std::vector<double> extractCurvatureOnPoints(const VCM& vcm);

	template <typename BackgroundPredicate, typename DTL2>
	DGtal::Z3i::DigitalSet extractSaddlePoints(const DGtal::Z3i::DigitalSet& setVolume,
											   const DTL2& dt,
											   const BackgroundPredicate& backgroundPredicate,
											   double R,
											   double r,
											   double delta);
};

template <typename VCM, typename Domain>
DGtal::Z3i::DigitalSet SaddleComputer::computeBranchingPartsWithVCMFeature(const VCM& vcm,
																	const Domain& domain, double threshold) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	DGtal::Z3i::DigitalSet aSet(domain);
	double R = vcm.R();
	for (P2EConstIterator  it = vcm.mapPoint2ChiVCM().begin(),
			  itE = vcm.mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
		if (ratio > threshold &&  sqrt(lambda[2]) < R*R)
			aSet.insert(it->first);
	}
	return aSet; 
}

template <typename VCM>
std::vector<double> SaddleComputer::extractCurvatureOnPoints(const VCM& vcm) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	std::vector<double> curvatureValues;
	for (P2EConstIterator  it = vcm.mapPoint2ChiVCM().begin(),
			  itE = vcm.mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
	    curvatureValues.push_back(ratio);
	}
	return curvatureValues; 
}

template <typename BackgroundPredicate, typename DTL2>
DGtal::Z3i::DigitalSet SaddleComputer::extractSaddlePoints(const DGtal::Z3i::DigitalSet& setVolume,
									const DTL2& dt,
									const BackgroundPredicate& backgroundPredicate,
									double R,
									double r,
									double delta) {
	typedef DGtal::ExactPredicateLpSeparableMetric<DGtal::Z3i::Space, 2> Metric;	
	typedef DGtal::functors::BallConstantPointFunction<DGtal::Z3i::Point, double> KernelFunction;

	typedef DGtal::Z3i::KSpace KSpace;
	typedef DGtal::ImplicitDigitalSurface< KSpace, BackgroundPredicate > DigitalSurfaceContainer;

	typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric, KernelFunction, DTL2> VCMOnSurface;
	typedef KSpace::Surfel Surfel;
	
	auto domainVolume = setVolume.domain();
	Metric l2;
	KSpace ks;
	ks.init( domainVolume.lowerBound(),
			 domainVolume.upperBound(), true );
	SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
	Surfel bel = Surfaces<KSpace>::findABel( ks, backgroundPredicate, 1000000 );
	DigitalSurfaceContainer* container =
		new DigitalSurfaceContainer( ks, backgroundPredicate, surfAdj, bel, false  );
	DigitalSurface< DigitalSurfaceContainer > surface( container ); //acquired

	//! [DVCM3D-instantiation]
	Surfel2PointEmbedding embType = InnerSpel; // Could be Pointels|InnerSpel|OuterSpel;
	KernelFunction chiSurface( 1.0, r );             // hat function with support of radius r

	
	VCMOnSurface* vcm_surface = new VCMOnSurface( surface, embType, R, r,
												  chiSurface, dt, delta, l2, true);
	
	
	std::vector<double> curvaturePoints = extractCurvatureOnPoints(*vcm_surface);
	double threshold2 = Statistics::unimodalThresholding(curvaturePoints);

//	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);
	Z3i::DigitalSet setSurface(domainVolume);
	std::set<Z3i::Point> pointSet;
	for ( auto it = surface.begin(), itE = surface.end(); it != itE; ++it )
		vcm_surface->getPoints( std::inserter( pointSet, pointSet.begin() ), *it );
	setSurface.insert(pointSet.begin(), pointSet.end());
	Z3i::DigitalSet branchingPoints = computeBranchingPartsWithVCMFeature(*vcm_surface, domainVolume, threshold2);
	return branchingPoints;
}

#endif
