#ifndef SLICE_UTILS_H
#define SLICE_UTILS_H

#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/images/ImageHelper.h"
//! [ExampleViewer3D2DImagesExtractImagesNonSliceHeader]
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/images/IntervalForegroundPredicate.h"
#include "Distance.h"

using namespace DGtal;

namespace SliceUtils {
	double deduceCoordinateFromEquationThreeUnknown(double a, double b, double c, double first_value, double second_value) {
		return (-(a * first_value + b * second_value) / c); 
	}

	template <typename Vector, typename Point>
	std::vector<Vector> computePlaneFromNormalVector(Vector normal, Point origin);

	template <typename Vector>
	Vector computeProjectionOnSet(const Z3i::DigitalSet& aSet, const Z3i::Point& aPoint,
								  const Z3i::RealPoint& aNormal);
	
	
	template <typename Pencil, typename Image>
	void slicesFromPlanes(Viewer3D<>&, const std::vector<Pencil> &, const Image&, std::string);


	template <typename ImageAdapter>
	double computeRadiusFromImage(const ImageAdapter& image, int thresholdMin, int thresholdMax);

}


template <typename Vector, typename Point>
std::vector<Vector> SliceUtils::computePlaneFromNormalVector(Vector normal, Point origin) {
	double a = normal[0];
	double b = normal[1];
	double c = normal[2];
	double d = -a*origin[0] - b*origin[1] - c*origin[2];

	Vector pRefOrigin;
	if ( a != 0 ) {
        pRefOrigin [0]= -d/a;
        pRefOrigin [1]= 0.0;
        pRefOrigin [2]= 0.0;
        if(pRefOrigin == origin) {
			pRefOrigin[1]=-1.0;
        }
	}
	else if ( b != 0 ) {
        pRefOrigin [0]= 0.0;
        pRefOrigin [1]= -d/b;
        pRefOrigin [2]= 0.0;
        if(pRefOrigin == origin) {
			pRefOrigin[0]=-1.0;
        }
	}
	else if ( c != 0 ) {
        pRefOrigin [0]= 0.0;
        pRefOrigin [1]= 0.0;
        pRefOrigin [2]= -d/c;
        if(pRefOrigin == origin) {
			pRefOrigin[0]=-1.0;
        }
	}
    Vector uDir1;
	uDir1=(pRefOrigin - origin)/((pRefOrigin - origin).norm());
	Vector uDir2;
	uDir2[0] = uDir1[1]*c-uDir1[2]*b;
	uDir2[1] = uDir1[2]*a-uDir1[0]*c;
	uDir2[2] = uDir1[0]*b-uDir1[1]*a;
      
	uDir2/=uDir2.norm();

	
	Vector myFirstAxisEmbeddedDirection = -uDir1;
	Vector mySecondAxisEmbeddedDirection = -uDir2;
	
	std::vector<Vector> fourPointsForPlane;
	Vector p1, p2, p3, p4;
	p1 = {(Vector)origin - myFirstAxisEmbeddedDirection - mySecondAxisEmbeddedDirection};
	p2 = {(Vector)origin - myFirstAxisEmbeddedDirection + mySecondAxisEmbeddedDirection};
	p3 = {(Vector)origin + myFirstAxisEmbeddedDirection + mySecondAxisEmbeddedDirection};
	p4 = {(Vector)origin + myFirstAxisEmbeddedDirection - mySecondAxisEmbeddedDirection};
	fourPointsForPlane = {p1, p2, p3, p4};
	return fourPointsForPlane;
}



/**
 * Constructs a distance transform from a 2D image
 * Returns the maximum distance of this DT in the object
 */
template <typename ImageAdapter>
double SliceUtils::computeRadiusFromImage(const ImageAdapter& image, int thresholdMin, int thresholdMax) {
    using namespace Z2i;

	typedef DGtal::functors::IntervalForegroundPredicate<ImageAdapter> Binarizer; 
	Binarizer binarizer(image, thresholdMin-1, thresholdMax);
	typedef DGtal::DistanceTransformation<Space, Binarizer, L2Metric> DTL2;
	
	DTL2 dt(&image.domain(), &binarizer, &l2Metric );
	double maxDT = (*std::max_element(dt.constRange().begin(), 
									  dt.constRange().end()));
	return maxDT;	
}

/**
 * Utility function designed at viewing and saving the 2D extracted slices
 */
template <typename Pencil, typename Image>
void SliceUtils::slicesFromPlanes(Viewer3D<>& viewer, const std::vector<Pencil> & vectorPlanes, const Image& volume, std::string outFileName) {
	typedef Image Image3D;
	typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D;
	typedef typename Image3D::Value Value;
	//! [ExampleViewer3D2DImagesExtractImagesNonSliceType]
	typedef DGtal::ConstImageAdapter<Image3D,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,Value, DGtal::functors::Identity> ImageAdapterExtractor;
	
	const int IMAGE_PATCH_WIDTH = 100;
	// Setting the image domain of the resulting image to be displayed in 3D:
	//! [ExampleViewer3D2DImagesExtractImagesNonSliceParam]
	unsigned int sliceNumber = 0;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH), volume.domain().upperBound() + Z3i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
									  DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::functors::Identity idV;
	
	for (auto it = ++vectorPlanes.rbegin(), itE = vectorPlanes.rend(); it != itE; ++it) {

		typename Pencil::Vector3d planeNormal = it->getTangent();
		typename Pencil::P origin = it->getPoint();
		DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, origin, planeNormal, IMAGE_PATCH_WIDTH, domain3Dyup.lowerBound());
		ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
		extractedImage.setDefaultValue(0);
		std::string outName;
		outName += outFileName + "_" + std::to_string(sliceNumber) + ".pgm";
		PGMWriter<ImageAdapterExtractor>::exportPGM(outName, extractedImage);
		viewer << extractedImage;
		viewer << DGtal::UpdateImage3DEmbedding<Z3i::Space, Z3i::KSpace>(sliceNumber, embedder(Z2i::RealPoint(0,0)), embedder(Z2i::RealPoint(IMAGE_PATCH_WIDTH,0)), embedder(domainImage2D.upperBound()), embedder(Z2i::RealPoint(0, IMAGE_PATCH_WIDTH)));
		sliceNumber++;
	}
}


#endif
