#ifndef SADDLE_COMPUTER_H
#define SADDLE_COMPUTER_H

#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "geometry/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "geometry/VCMUtil.h"


template <typename DTL2, typename BackgroundPredicate>
class SaddleComputer {
public:
	typedef DGtal::ExactPredicateLpSeparableMetric<DGtal::Z3i::Space, 2> Metric;
	typedef DGtal::functors::BallConstantPointFunction<DGtal::Z3i::Point, double> KernelFunction;

	typedef DGtal::Z3i::KSpace KSpace;
	typedef DGtal::ImplicitDigitalSurface< KSpace, BackgroundPredicate > DigitalSurfaceContainer;

	typedef KSpace::Surfel Surfel;
	typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric, KernelFunction, DTL2> VCMOnSurface;
	typedef DGtal::Z3i::Domain Domain;

public:

	SaddleComputer(const DGtal::Z3i::DigitalSet& setVolume,
				   const DTL2& dt,
				   const BackgroundPredicate& backgroundPredicate,
				   double R,
				   double r,
				   double delta);

	~SaddleComputer() {
		delete myVCMSurface;
		delete mySurface;
	}

	template <typename VCM>
	DGtal::Z3i::DigitalSet computeBranchingPartsWithVCMFeature(const VCM& vcm, double threshold);

	std::vector<double> extractCurvatureOnPoints();

	std::vector<std::pair<Z3i::Point, double> > mapPointToEigenvalue();

	DGtal::Z3i::DigitalSet extractSaddlePoints(const DGtal::Z3i::DigitalSet& setVolume);

	std::vector<DGtal::Z3i::Object26_6> saddleConnectedComponents(const DGtal::Z3i::DigitalSet& saddles);

	template <typename Matrix>
	DGtal::Z3i::DigitalSet saddlePointsToOnePoint(const std::vector<DGtal::Z3i::Object26_6>& branchingPoints);



private:
	VCMOnSurface* myVCMSurface;
	DigitalSurface< DigitalSurfaceContainer >* mySurface;
	Domain myDomain;
};

template <typename DTL2, typename BackgroundPredicate>
SaddleComputer<DTL2, BackgroundPredicate>::SaddleComputer(const DGtal::Z3i::DigitalSet& setVolume,
									const DTL2& dt,
									const BackgroundPredicate& backgroundPredicate,
									double R,
									double r,
									double delta) {


	myDomain = setVolume.domain();
	Metric l2;
	KSpace ks;
	ks.init( myDomain.lowerBound(),
			 myDomain.upperBound(), true );
	SurfelAdjacency<KSpace::dimension> surfAdj( true ); // interior in all directions.
	Surfel bel = Surfaces<KSpace>::findABel( ks, backgroundPredicate, 1000000 );
	DigitalSurfaceContainer* container =
		new DigitalSurfaceContainer( ks, backgroundPredicate, surfAdj, bel, false  );
    mySurface = new DigitalSurface< DigitalSurfaceContainer >( container ); //acquired

	//! [DVCM3D-instantiation]
	Surfel2PointEmbedding embType = Pointels; // Could be Pointels|InnerSpel|OuterSpel;
	KernelFunction chiSurface( 1.0, r );             // hat function with support of radius r


    myVCMSurface= new VCMOnSurface( *mySurface, embType, R, r,
									chiSurface, dt, delta, l2, true);
}

template <typename DTL2, typename BackgroundPredicate>
template <typename VCM>
DGtal::Z3i::DigitalSet SaddleComputer<DTL2, BackgroundPredicate>::computeBranchingPartsWithVCMFeature(const VCM& vcm, double threshold) {
	typedef typename VCM::Point2EigenStructure::const_iterator P2EConstIterator;

	DGtal::Z3i::DigitalSet aSet(myDomain);
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

template <typename DTL2, typename BackgroundPredicate>
std::vector<double> SaddleComputer<DTL2, BackgroundPredicate>::extractCurvatureOnPoints() {
	typedef typename VCMOnSurface::Point2EigenStructure::const_iterator P2EConstIterator;

	std::vector<double> curvatureValues;
	for (P2EConstIterator  it = myVCMSurface->mapPoint2ChiVCM().begin(),
			  itE = myVCMSurface->mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
	    curvatureValues.push_back(ratio);
	}
	return curvatureValues;
}


template <typename DTL2, typename BackgroundPredicate>
std::vector<std::pair<Z3i::Point, double> > SaddleComputer<DTL2, BackgroundPredicate>::mapPointToEigenvalue() {
	typedef typename VCMOnSurface::Point2EigenStructure::const_iterator P2EConstIterator;

	std::vector<std::pair<Z3i::Point, double> > curvatureValues;
	for (P2EConstIterator  it = myVCMSurface->mapPoint2ChiVCM().begin(),
			  itE = myVCMSurface->mapPoint2ChiVCM().end(); it != itE; ++it )
    {
		auto lambda = it->second.values;
		double ratio = VCMUtil::computeCurvatureJunction(lambda);
	    curvatureValues.push_back(std::make_pair(it->first, ratio));
	}
	return curvatureValues;
}

template <typename DTL2, typename BackgroundPredicate>
DGtal::Z3i::DigitalSet SaddleComputer<DTL2, BackgroundPredicate>::extractSaddlePoints(const DGtal::Z3i::DigitalSet& setVolume){

	auto myDomain = setVolume.domain();
	std::vector<double> curvaturePoints = extractCurvatureOnPoints();
	double threshold2 = Statistics::unimodalThresholding(curvaturePoints);

//	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);
	Z3i::DigitalSet branchingPoints = computeBranchingPartsWithVCMFeature(*myVCMSurface, threshold2);
	return branchingPoints;
}


/*
 * Reduce each set of saddle points computed from extractSaddlePoints to one point
 * Each point have the maximal curvature value within its connected component
 */
template <typename DTL2, typename BackgroundPredicate>
std::vector<DGtal::Z3i::Object26_6> SaddleComputer<DTL2, BackgroundPredicate>::saddleConnectedComponents(const DGtal::Z3i::DigitalSet& branchingPoints) {
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;

	auto myDomain = branchingPoints.domain();
	NotPointPredicate notBranching(branchingPoints);
	Z3i::Object26_6 obj(Z3i::dt26_6, branchingPoints);
	std::vector<DGtal::Z3i::Object26_6> objects;
	std::back_insert_iterator< std::vector<Z3i::Object26_6> > inserter( objects );
	unsigned int nbConnectedComponents = obj.writeComponents(inserter);
	return objects;

}

/*
 * Reduce each set of saddle points computed from extractSaddlePoints to one point
 * Each point have the maximal curvature value within its connected component
 */
template <typename DTL2, typename BackgroundPredicate>
template <typename Matrix>
DGtal::Z3i::DigitalSet SaddleComputer<DTL2, BackgroundPredicate>::saddlePointsToOnePoint(const std::vector<DGtal::Z3i::Object26_6>& objects) {
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;

	Z3i::DigitalSet maxCurvaturePoints(myDomain);
	Matrix vcmrB, evecB;
	Z3i::RealVector evalB;
	for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
		double ratioMax = 0;
		Z3i::Point maximizingCurvaturePoint;
		for (auto itPoint = it->pointSet().begin(), itPointE = it->pointSet().end(); itPoint != itPointE; ++itPoint) {
			auto lambda = (myVCMSurface->mapPoint2ChiVCM()).at(*itPoint).values;
			double ratio = VCMUtil::computeCurvatureJunction(lambda);
			if (ratio > ratioMax) {
				ratioMax = ratio;
				maximizingCurvaturePoint = *itPoint;
			}
		}
		maxCurvaturePoints.insert(maximizingCurvaturePoint);
	}
	return maxCurvaturePoints;
}



#endif
