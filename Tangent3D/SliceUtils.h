#ifndef SLICE_UTILS_H
#define SLICE_UTILS_H

#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageHelper.h"
//! [ExampleViewer3D2DImagesExtractImagesNonSliceHeader]
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "../CreateCurve/Distance.h"

using namespace DGtal;

namespace SliceUtils {
	double deduceCoordinateFromEquationThreeUnknown(double a, double b, double c, double first_value, double second_value) {
		return (-(a * first_value + b * second_value) / c); 
	}

	template <typename Vector>
	std::vector<Vector> computePlaneFromNormalVector(Vector normal);
	
	template <typename Pencil, typename Image>
	void slicesFromPlanes(Viewer3D<>&, const std::vector<Pencil> &, const Image&, std::string);
}

template <typename Vector>
std::vector<Vector> SliceUtils::computePlaneFromNormalVector(Vector normal) {
	double a = normal[0];
	double b = normal[1];
	double c = normal[2];
	std::vector<Vector> fourPointsForPlane;

	double x, y, z;
	Vector p1, p2, p3, p4;
	if (a != 0) {
		y = 1;
		z = 1;
		x = SliceUtils::deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p1 = Vector(x, y, z).getNormalized();
		y = -1;
		x = SliceUtils::deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p2 = Vector(x, y, z).getNormalized();
		z = -1;
		x = SliceUtils::deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p3 = Vector(x, y, z).getNormalized();
		y = 1;
		x = SliceUtils::deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p4 = Vector(x, y, z).getNormalized();
	} else if (b != 0) {
		x = 1;
		z = 1;
		y = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p1 = Vector(x, y, z).getNormalized();
		x = -1;
		y = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p2 = Vector(x, y, z).getNormalized();
		z = -1;
		y = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p3 = Vector(x, y, z).getNormalized();
		x = 1;
		y = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p4 = Vector(x, y, z).getNormalized();
	} else if (c != 0) {
		x = 1;
		y = 1;
		z = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p1 = Vector(x, y, z).getNormalized();
		y = -1;
		z = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p2 = Vector(x, y, z).getNormalized();
		x = -1;
		z = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p3 = Vector(x, y, z).getNormalized();
		y = 1;
		z = SliceUtils::deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p4 = Vector(x, y, z).getNormalized();
	}
	fourPointsForPlane = {p1, p2, p3, p4};
    
	return fourPointsForPlane;
}

template <typename Pencil, typename Image>
void SliceUtils::slicesFromPlanes(Viewer3D<>& viewer, const std::vector<Pencil> & vectorPlanes, const Image& volume, std::string outFileName) {
	typedef Image Image3D;
	typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> Image2D;
	typedef typename Image3D::Value Value;
	//! [ExampleViewer3D2DImagesExtractImagesNonSliceType]
	typedef DGtal::ConstImageAdapter<Image3D,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,Value, DGtal::functors::Identity> ImageAdapterExtractor;
	
	const int IMAGE_PATCH_WIDTH = 100;
	// Setting the image domain of the resulting image to be displayed in 3D:
	//DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
	//								  DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH)); 
	//! [ExampleViewer3D2DImagesExtractImagesNonSliceParam]
	unsigned int sliceNumber = 0;
	Z3i::Domain domain3Dyup(volume.domain().lowerBound() + Z3i::Point(-IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH, -IMAGE_PATCH_WIDTH), volume.domain().upperBound() + Z3i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
											  DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH));
	DGtal::functors::Identity idV;
	for (auto it = vectorPlanes.begin(), itE = vectorPlanes.end(); it != itE; ++it) {

		typename Pencil::Vector3d planeNormal = it->getTangent();
		typename Pencil::P origin = it->getPoint();
		std::vector<typename Pencil::Vector3d> pointsForPlane = SliceUtils::computePlaneFromNormalVector(planeNormal);	 
		DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, origin,planeNormal, IMAGE_PATCH_WIDTH);
		ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
		// Image2D image(domainImage2D);
		// Z2i::Point center = {50, 50};
		// for (auto it = extractedImage.domain().begin(), itE = extractedImage.domain().end(); it != itE; ++it) {
		// 	if (euclideanDistance(*it, center) < 20.0 && extractedImage(*it) > 0)
		// 	{
		// 		image.setValue(*it, 1);
		// 	}
		// }
		std::string outName;
		outName += outFileName + "_" + std::to_string(sliceNumber) + ".pgm";
		GenericWriter<ImageAdapterExtractor>::exportFile(outName, extractedImage);
		if (sliceNumber < 88) {
			viewer << extractedImage;
			viewer << DGtal::UpdateImage3DEmbedding<Z3i::Space, Z3i::KSpace>(sliceNumber, embedder(Z2i::RealPoint(0,0)), embedder(Z2i::RealPoint(IMAGE_PATCH_WIDTH,0)), embedder(domainImage2D.upperBound()), embedder(Z2i::RealPoint(0, IMAGE_PATCH_WIDTH)));
		}
		sliceNumber++;
	}
}

#endif
