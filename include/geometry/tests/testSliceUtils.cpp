#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/PGMWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "../ImageUtil.h"
#include "../SliceUtils.h"
#include "../../Statistics.h"
#include "../PointUtil.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "../../shapes/Ball.h"
#include "../VCMUtil.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"


using namespace DGtal;
using namespace std;
Z3i::DigitalSet computeHalfShell(const Z3i::Point& center, const Z3i::RealVector& dirVector,
								 const Z3i::DigitalSet& plane,
								 const Z3i::DigitalSet& setVolume, double radiusInnerBall) {
	typedef BreadthFirstVisitor<Z3i::Object26_6, std::set<Z3i::Point> > Visitor;
	typedef Visitor::Node Node;

	Z3i::Object26_6 obj(Z3i::dt26_6, setVolume);
    Visitor visitor(obj, plane.begin(), plane.end());
	Z3i::DigitalSet shell(setVolume.domain());
	int radius = (int) radiusInnerBall;
	while (!visitor.finished()) {
		Node node = visitor.current();
		if (node.second > radius) break;
		if (node.second == radius && VCMUtil::abovePlane(node.first, dirVector, center))
			shell.insert(node.first);
		visitor.expand();

	}
    return shell;
}

unsigned int computeDegree(const Z3i::DigitalSet& shell) {
	typedef Z3i::Object26_6 ObjectType;

	ObjectType objectImage(Z3i::dt26_6, shell);
	vector<ObjectType> objects;
	back_insert_iterator< std::vector<ObjectType> > inserter( objects );
	unsigned int nbConnectedComponents = objectImage.writeComponents(inserter);
	unsigned int cpt = 0;

	for (const auto& obj : objects) {
		if (obj.size() > 1)
			cpt++;
	}
	return cpt;
}

void testSliceFromPlane(int argc, char** argv) {

	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;
	typedef ImageSelector<Z2i::Domain, unsigned char>::Type Image2D;

	QApplication app(argc, argv);
	Viewer3D<> viewer;
	viewer.show();

	Image image = GenericReader<Image>::import("/home/florent/test_img/bronche.vol");
	Z3i::DigitalSet setVolume(image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image>(setVolume, image, 0, 255);
	Z3i::Point origin(204,235,329);
	Z3i::RealPoint normal(-0.740353,0.0560496,-0.66987);
    Z2i::Point center(50,50);
	ImageAdapterExtractor extractedImage = SliceUtils::sliceFromPlane<ImageAdapterExtractor>(-normal, origin, image, 100);
	Image2D processImage = ImageUtil::convertImage<Image2D>(extractedImage);
	auto domainImage2D = processImage.domain();
	Z2i::DigitalSet aSet(domainImage2D);
	for (const Z2i::Point& p : domainImage2D) {
		if (processImage(p) >= 1) {
			aSet.insert(p);
		}
	}
	//	PGMWriter<ImageAdapterExtractor>::exportPGM("slice.pgm", extractedImage);
    Eigen::MatrixXd covmatrix = Statistics::computeCovarianceMatrix<Eigen::MatrixXd>(aSet);
	if (covmatrix.size() == 0) return;
    Z2i::RealVector projection = Statistics::extractEigenVector<Z2i::RealVector>(covmatrix, 1);
    Z2i::Point trackedPoint = PointUtil::trackPoint(center, aSet, projection);
    Z2i::Point otherTrackedPoint = PointUtil::trackPoint(center, aSet, -projection);
	double distance = Z2i::l2Metric(otherTrackedPoint, trackedPoint) / 2.0;
    DGtal::trace.info() << distance << endl;

	Ball<Z3i::Point> ball(origin, distance);
	Z3i::DigitalSet startingPoint(setVolume.domain());
	startingPoint.insert(origin);
	Z3i::DigitalSet pointInBall = computeHalfShell(origin, Z3i::RealVector(0,0,0), startingPoint, setVolume, distance);
	unsigned int nbCC = computeDegree(pointInBall);
	DGtal::trace.info() << nbCC << endl;
	for (const Z3i::Point& p : pointInBall)
		viewer << CustomColors3D(Color::Green, Color::Green) << p;
	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << origin;
	viewer << CustomColors3D(Color(128,128,128,20), Color(128,128,128,20)) << setVolume;
	viewer << Viewer3D<>::updateDisplay;


	app.exec();

}


int main(int argc, char** argv) {
	testSliceFromPlane(argc, argv);
	return 0;
}
