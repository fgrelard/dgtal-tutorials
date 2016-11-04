#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "geometry/PointUtil.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
#include "DGtal/io/viewers/Viewer3D.h"

using namespace DGtal;
using namespace std;

void testBezierDeCasteljau(Viewer3D<>& viewer) {
	Z3i::Point first(2, 16, 31);
	Z3i::Point destination(3, -1, 23);
	Z3i::Point control1(0, 1, 25);
	Z3i::Point control2(1, 5, 25);

	vector<Z3i::Point> points = PointUtil::bezierCurveDeCasteljau(first, destination, control1, control2);
	for (const Z3i::Point& p : points) {
		viewer << CustomColors3D(Color::Red, Color::Red) << p;
	}
	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << first << destination << control1 << control2;
	viewer << Viewer3D<>::updateDisplay;
}

void testLinking() {
	using namespace Z3i;
	Point first(233, 276, 172);
	Point second(226, 276, 185);

	vector<Point> points = PointUtil::linkTwoPoints(first, second);
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		trace.info() << *it << endl;
	}
}

void testLinking26(Viewer3D<>& viewer) {
	using namespace Z3i;
	Point first(0, 0, 0);
	Point second(-1, -3, -9);

	vector<Point> points = PointUtil::linkTwoPoints26(first, second);
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		trace.info() << *it << endl;
	}
	for (const Z3i::Point& p : points) {
		viewer << CustomColors3D(Color::Red, Color::Red) << p;
	}
	viewer << CustomColors3D(Color::Yellow, Color::Yellow) << first << second;
	viewer << Viewer3D<>::updateDisplay;
}

void testDSSLinking() {
	typedef vector<Z3i::Point>::iterator Iterator;
	typedef StandardDSS6Computer<Iterator, int, 8> DSS;
	typedef SaturatedSegmentation<DSS> Segmentation;

	vector<Z3i::Point> vPoints;
	vPoints.push_back(Z3i::Point(233, 276, 172));
	vPoints.push_back(Z3i::Point(226, 276, 185));

	DSS algo;
	Iterator i = vPoints.begin();
	algo.init(i);
	trace.info() << "init with " << (*i) << std::endl;

    while (algo.extendFront()) {
      trace.info() << "extended with " << (*(--algo.end())) << std::endl;
    }

}

int main(int argc, char** argv) {
	QApplication app(argc, argv);
	Viewer3D<> viewer;
	viewer.show();
	testLinking26(viewer);
//	testBezierDeCasteljau(viewer);
//	testLinking();
//	testDSSLinking();
	app.exec();

	return 0;
}
