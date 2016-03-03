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

using namespace DGtal;
using namespace std;

void testSurfaceTraversal(int argc, char** argv) {
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef Z3i::Object26_6 Graph;
	typedef DistanceToPointFunctor<Z3i::L2Metric> Distance;
	typedef DistanceBreadthFirstVisitor<Graph, Distance, std::set<Z3i::Point>> Visitor;
	
	Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");
	Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(image, 1, 255);
	Graph graphSurface(Z3i::dt26_6, setSurface);

	Z3i::Point source(1, -9, 2);
	Z3i::Point destination(-9, -1, 8);
	vector<Z3i::Point> path = SurfaceTraversal::AStarAlgorithm(graphSurface, source, destination);
	//vector<Z3i::Point> path = SurfaceTraversal::DijkstraAlgorithm(graphSurface, source, destination);
	
	QApplication app(argc, argv);
	Viewer3D<> viewer;
	viewer.show();
	const Color CURVE3D_COLOR( 100, 100, 140, 128 );
   
	for (auto it = path.begin(), ite = path.end(); it != ite; ++it) {
		viewer << CustomColors3D(Color::Green, Color::Green) << *it;
	}
   	
	viewer << CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR) << graphSurface;
	viewer << Viewer3D<>::updateDisplay;
	app.exec();
}
	
int main(int argc, char** argv) {
	testSurfaceTraversal(argc, argv);
	return 0;
}
