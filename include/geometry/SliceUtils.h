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
#include <Eigen/Dense>

using namespace DGtal;

namespace SliceUtils {
	double deduceCoordinateFromEquationThreeUnknown(double a, double b, double c, double first_value, double second_value) {
		return (-(a * first_value + b * second_value) / c); 
	}

	template <typename Vector, typename Point>
	std::vector<Vector> computePlaneFromNormalVector(Vector normal, Point origin);

	template <typename Vector, typename Point>
	Vector computeNormalFromLinearRegression(const std::vector<Point>& points);

	template <typename Vector, typename Point>
	Vector computeNormalFromCovarianceMatrix(const std::vector<Point>& points);

	template <typename Image2D>
	Z2i::RealPoint centerOfMass(const Image2D& image);
	
	template <typename Pencil, typename Image>
	void slicesFromPlanes(Viewer3D<>&, const std::vector<Pencil> &, const Image&, std::string);


	/**
	 * DO NOT USE
	 * Does not work, doesnt take into account the actual pixel values
	 */
	template <typename ImageAdapter, typename Vector, typename Point, typename Image>
	ImageAdapter sliceFromPlane(const Vector& normal, const Point& origin, const Image& image, const int patch_width);
	
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
 * Computes the normal of a plane from a set of points
 * Method : linear regression
 */
template <typename Vector, typename Point>
Vector SliceUtils::computeNormalFromLinearRegression(const std::vector<Point> & points) {
	typedef Eigen::Matrix<double, Eigen::Dynamic, 3> MatrixXi;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
	unsigned int size = points.size();
	MatrixXi A(size, 3);
	VectorXi b = VectorXi::Zero(size, 1);
	
	for (int i = 0; i < size; i++) {
		A(i, 0) = (double)points[i][0]*1.0;
		A(i, 1) = (double)points[i][1]*1.0;
		A(i, 2) = 1.0;
		b(i, 0) = (double)points[i][2]*1.0;
	}
	Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
	Vector normal;
	normal[0] = x(0, 0);
	normal[1] = x(1, 0);
	normal[2] = -1;
	return normal.getNormalized();
}

/**
 * Computes the normal of a plane from a set of points
 * Method : covariance matrix
 */
template <typename Vector, typename Point>
Vector SliceUtils::computeNormalFromCovarianceMatrix(const std::vector<Point> & points) {
	typedef Eigen::MatrixXd MatrixXd;
	
	unsigned int size = points.size();
	if (size < 2) return Vector::zero;
	
	MatrixXd A(size, 3);
	for (int i = 0; i < size; i++) {
		A(i, 0) = (double)points[i][0] * 1.0;
		A(i, 1) = (double)points[i][1] * 1.0;
		A(i, 2) = (double)points[i][2] * 1.0;
	}
	MatrixXd centered = A.rowwise() - A.colwise().mean();
	MatrixXd cov = (centered.adjoint() * centered) / double(A.rows() - 1);
	Eigen::SelfAdjointEigenSolver<MatrixXd> eig(cov);
	Vector normal;
	auto veigen = eig.eigenvectors().col(0);
	normal[0] = veigen[0];
	normal[1] = veigen[1];
	normal[2] = veigen[2];
	return normal;
}

/**
 * Constructs the intersection of a given object image with a plane (defined by normal and origin)
 * Returns the corresponding 2D intersection
 */
template <typename ImageAdapter, typename Vector, typename Point, typename Image>
ImageAdapter SliceUtils::sliceFromPlane(const Vector& normal, const Point& origin, const Image& image, const int patch_width) {
	const Z3i::Domain domain3Dyup(image.domain().lowerBound() + Z3i::Point(-patch_width, -patch_width, -patch_width), image.domain().upperBound() + Z3i::Point(patch_width, patch_width, patch_width));
	const DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
											DGtal::Z2i::Point(patch_width, patch_width));
	DGtal::functors::Identity idV;
	DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(domain3Dyup, origin, normal, patch_width);

	ImageAdapter extractedImage(image, domainImage2D, embedder, idV);
	return extractedImage;
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

template <typename Image2D>
Z2i::RealPoint SliceUtils::centerOfMass(const Image2D& image) {
	double m00 = 0;
	double m10 = 0;
	double m01 = 0;

    for (auto it = image.domain().begin(), ite = image.domain().end(); it != ite; ++it) {
		Z2i::Point current = *it;
		m00 += image(current);
		m10 += current[0] * image(current);
		m01 += current[1] * image(current);
	}
	if (m00 != 0) 
		return Z2i::RealPoint(m10/m00, m01/m00);
	else
		return Z2i::RealPoint();
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
