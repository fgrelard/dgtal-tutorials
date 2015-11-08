#include "../Statistics.h"
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/GenericReader.h"

using namespace DGtal;
using namespace std;

void testCenterOfMass() {
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;

	Image image = GenericReader<Image>::import("/home/florent/test_img/slices/boudin/slice_1.pgm");
	
	Z2i::Point point = Statistics::centerOfMass(image);
	trace.info() << point << endl;
}

void testComputeCovarianceMatrix() {
	typedef Eigen::MatrixXd MatrixXd;
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image;
	typedef DGtal::Z2i::RealPoint RealPoint;
	
	Image image = GenericReader<Image>::import("/home/florent/trash/slice_250.pgm");
	MatrixXd matrixCovariance = Statistics::computeCovarianceMatrix<MatrixXd>(image);
	RealPoint vectorZero = Statistics::extractEigenVector<RealPoint>(matrixCovariance, 0);
	RealPoint vectorOne = Statistics::extractEigenVector<RealPoint>(matrixCovariance, 1);
	trace.info() << vectorZero << " " << vectorOne << endl;
}

int main() {
	testCenterOfMass();
	testComputeCovarianceMatrix();
	return 0;
}
