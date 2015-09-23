#ifndef __STATISTICS__
#define __STATISTICS
#include <vector>
#include <algorithm>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include <Eigen/Dense>

namespace Statistics {
	template <typename T>
	double mean(const std::vector<T>& aVector);

	template <typename T>
	double stddev(const std::vector<T>& aVector);

	template <typename Image2D>
	DGtal::Z2i::RealPoint centerOfMass(const Image2D& image);

	template <typename Image3D>
	DGtal::Z3i::RealPoint centerOfMass3D(const Image3D& image);
	
	DGtal::Z2i::RealPoint extractCenterOfMass(const DGtal::Z2i::DigitalSet& set);
	
	DGtal::Z3i::RealPoint extractCenterOfMass3D(const DGtal::Z3i::DigitalSet& set);
	
	template <typename Vector, typename Point>
	Vector computeNormalFromLinearRegression(const std::vector<Point>& points);

	template <typename Vector, typename Point>
	Vector computeNormalFromCovarianceMatrix(const std::vector<Point>& points);

	template <typename Vector, typename Matrix>
	Vector extractEigenVector(const Matrix& m, int colNumber);

	template <typename Matrix, typename Image2D>
	Matrix computeCovarianceMatrix(const Image2D& image);
	
}

template <typename T>
double Statistics::mean(const std::vector<T>& aVector) {
	return std::accumulate( aVector.begin(), aVector.end(), 0.0f )/ aVector.size();
}

template <typename T>
double Statistics::stddev(const std::vector<T>& aVector) {
	std::vector<double> zero_mean( aVector );
	transform( zero_mean.begin(), zero_mean.end(), zero_mean.begin(),bind2nd( std::minus<double>(), mean(aVector) ) );

	double deviation = inner_product( zero_mean.begin(),zero_mean.end(), zero_mean.begin(), 0.0f );
	deviation = sqrt( deviation / ( aVector.size() - 1 ) );
	return deviation;
}

template <typename Image2D>
DGtal::Z2i::RealPoint Statistics::centerOfMass(const Image2D& image) {
	double m00 = 0;
	double m10 = 0;
	double m01 = 0;

	for (auto it = image.domain().begin(), ite = image.domain().end(); it != ite; ++it) {
		DGtal::Z2i::Point current = *it;
		m00 += image(current);
		m10 += current[0] * image(current);
		m01 += current[1] * image(current);
	}
	if (m00 != 0) 
		return DGtal::Z2i::RealPoint(m10/m00, m01/m00);
	else
		return DGtal::Z2i::RealPoint();
}

template <typename Image3D>
DGtal::Z3i::RealPoint Statistics::centerOfMass3D(const Image3D& image) {
	double m000 = 0;
	double m100 = 0;
	double m010 = 0;
	double m001 = 0;

	for (auto it = image.domain().begin(), ite = image.domain().end(); it != ite; ++it) {
		DGtal::Z3i::Point current = *it;
		m000 += image(current);
		m100 += current[0] * image(current);
		m010 += current[1] * image(current);
		m001 += current[2] * image(current);
	}
	if (m000 != 0) 
		return DGtal::Z3i::RealPoint(m100/m000, m010/m000, m001/m000);
	else
		return DGtal::Z3i::RealPoint();
}

	
DGtal::Z2i::RealPoint Statistics::extractCenterOfMass(const DGtal::Z2i::DigitalSet& set) {
	if (set.size() != 0) {
		typedef DGtal::ImageSelector<DGtal::Z2i::Domain, unsigned char>::Type Image2D;
		Image2D image2D = DGtal::ImageFromSet<Image2D>::create(set, 150);
		DGtal::Z2i::RealPoint centerOfMass = Statistics::centerOfMass(image2D);
		return centerOfMass;
	}
	else
		return DGtal::Z2i::RealPoint();
}

DGtal::Z3i::RealPoint Statistics::extractCenterOfMass3D(const DGtal::Z3i::DigitalSet& set) {
	if (set.size() != 0) {
		typedef DGtal::ImageSelector<DGtal::Z3i::Domain, unsigned char>::Type Image3D;
		Image3D image3D = DGtal::ImageFromSet<Image3D>::create(set, 150);
		DGtal::Z3i::RealPoint centerOfMass = centerOfMass3D(image3D);
		return centerOfMass;
	}
	else
		return DGtal::Z3i::RealPoint();
}

/**
 * Computes the normal of a plane from a set of points
 * Method : linear regression
 */
template <typename Vector, typename Point>
Vector Statistics::computeNormalFromLinearRegression(const std::vector<Point> & points) {
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
Vector Statistics::computeNormalFromCovarianceMatrix(const std::vector<Point> & points) {
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

template <typename Matrix, typename Image2D>
Matrix Statistics::computeCovarianceMatrix(const Image2D& image) {
	typedef typename Image2D::Domain Domain;
	typedef typename Domain::Point Point;
	int size = 0;
	for (typename Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
		 it != ite; ++it) {
		Point point = *it;
		if (image(*it) > 0) {
			size++;
		}
	}
	
	Matrix A(size, 2);
	int i = 0;
	for (typename Domain::ConstIterator it = image.domain().begin(), ite = image.domain().end();
		 it != ite; ++it) {
		Point point = *it;
		if (image(*it) > 0) {
			A(i, 0) = (double) point[0] * 1.0;
			A(i, 1) = (double) point[1] * 1.0;
			i++;
		}
	}
	Matrix centered = A.rowwise() - A.colwise().mean();
	Matrix cov = (centered.adjoint() * centered) / double(A.rows() - 1);
    return cov;
}

template <typename Vector, typename Matrix>
Vector Statistics::extractEigenVector(const Matrix& m, int colNumber) {
	Eigen::SelfAdjointEigenSolver<Matrix> eig(m);
	Vector normal;
	auto veigen = eig.eigenvectors().col(colNumber);
	normal[0] = veigen[0];
	normal[1] = veigen[1];
	normal[2] = veigen[2];
	return normal;
}


#endif
