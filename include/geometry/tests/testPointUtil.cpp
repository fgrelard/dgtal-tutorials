#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "geometry/PointUtil.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h" 
#include "DGtal/geometry/curves/ArithmeticalDSS.h"

using namespace DGtal;
using namespace std;

void testLinking() {
	using namespace Z3i;
	Point first(233, 276, 172);
	Point second(226, 276, 185);

	vector<Point> points = PointUtil::linkTwoPoints(first, second);
	for (auto it = points.begin(), ite = points.end(); it != ite; ++it) {
		trace.info() << *it << endl;
	}
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

int main() {
	testLinking();
//	testDSSLinking();
	return 0;
}
