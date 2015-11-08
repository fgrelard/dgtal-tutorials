#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/graph/DistanceBreadthFirstVisitor.h"
#include "DGtal/io/viewers/Viewer3D.h"

#include "surface/SurfaceUtils.h"
#include "surface/SurfaceTraversal.h"
#include "geometry/DistanceToPointFunctor.h"

#include <QtGui/qapplication.h>
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "geometry/GeodesicMetric.h"

using namespace DGtal;
using namespace std;

void testGeodesicMetric(int argc, char** argv) {
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef Z3i::Object26_6 Graph;
	typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
	typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;
	typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;
	typedef GeodesicMetric<Z3i::Space, Graph> GeodesicMetric;
	typedef Z3i::L2Metric Metric;
	typedef VoronoiMap<Z3i::Space, NotPointPredicate, Metric> VoronoiMap;
	Image image = GenericReader<Image>::import("/home/florent/test_img/bronche2.vol");
	Z3i::Domain domainVolume = image.domain();
	Z3i::DigitalSet setVolume(domainVolume);
	SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, image, 
												  0, 255);
	Graph graphSurface(Z3i::dt26_6, setVolume);

	Z3i::Point source(10, 0, 0);
	Z3i::DigitalSet aSet(domainVolume);
	aSet.insert(Z3i::Point(250,250,250));
	aSet.insert(Z3i::Point(350,200,100));
	NotPointPredicate notPredicate(aSet);
	Metric l2;
	GeodesicMetric geoMetric(graphSurface);
//	VoronoiMap vMap(domainVolume, notPredicate, geoMetric);
	VoronoiMap vMap(domainVolume, notPredicate, l2);
	QApplication app(argc, argv);
	Viewer3D<> viewer;
	viewer.show();
	const Color CURVE3D_COLOR( 100, 100, 140, 128 );


	for (auto it = setVolume.begin(), ite = setVolume.end(); it != ite; ++it) {
		if (vMap(*it) == Z3i::Point(250,250,250)) {
			viewer << CustomColors3D(Color::Red, Color::Red) << *it;
		}
		else {
			viewer << CustomColors3D(Color::Blue, Color::Blue) << *it;
		}
	}

   	
	viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << graphSurface;
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
}
	
int main(int argc, char** argv) {
	testGeodesicMetric(argc, argv);
	return 0;
}
