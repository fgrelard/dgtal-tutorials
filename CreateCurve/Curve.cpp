/** File allowing to create a curve in DGtal according to an equation **/
#include <iostream>
#include <fstream>
#include <time.h>
#include "DGtal/base/Common.h"
#include "DGtal/io/Display3DFactory.h"
// Shape construction
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
// Drawing
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
#include "DGtal/io/writers/VolWriter.h"
#include "Distance.h"
#include "Ball.h"
#include "ModelisationUtils.h"
///////////////////////////////////////////////////////////////////////////////
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/io/boards/Board2D.h"
#include "../Tangent3D/Pencil.h"
#include "../Tangent3D/MSTTangent.h"

using namespace std;
using namespace DGtal;

template <typename Point>
class PointDeletable : public Point {
public:
	typedef typename Point::Component Scalar;
public:
	PointDeletable() : Point() {}
	PointDeletable(const Scalar& x, const Scalar& y, const Scalar& z) : Point(x,y,z) {}
	bool myMarkedToDelete = false;
};


typedef Z3i::DigitalSet DigitalSet;
typedef Z3i::Space Space;
typedef Z3i::KSpace KSpace;
typedef Z3i::Curve::ArrowsRange Range;
typedef PointDeletable<Z3i::Point> Point;
typedef Range::ConstIterator ConstIterator;
typedef vector<Point>::const_iterator PointIterator;


void create2DNaive() {
	Z2i::Domain domain({0,0}, {10,10});
	Z2i::DigitalSet dSet(domain);
	Board2D board;
 	create2DCurve(dSet);
	//create2DNaiveTangentsForVisu(dSet, board);

	typedef Z2i::DigitalSet::Iterator Iterator;
	typedef typename IteratorCirculatorTraits<Iterator>::Value::Coordinate Coordinate; 
	typedef ArithmeticalDSSComputer<Iterator,Coordinate,8> RecognitionAlgorithm;
	typedef SaturatedSegmentation<RecognitionAlgorithm> Segmentation;
   
	RecognitionAlgorithm algo;
	Segmentation s(dSet.begin(),dSet.end(),algo);
     
    typename Segmentation::SegmentComputerIterator i = s.begin();
	typename Segmentation::SegmentComputerIterator end = s.end();
	board << domain;
	for (; i != end; ++i) {
         
		typename Segmentation::SegmentComputerIterator::SegmentComputer segment(*i); 
		
       board << SetMode(segment.primitive().className(), "Points") << segment.primitive() << SetMode(segment.primitive().className(), "BoundingBox" )
			 << segment.primitive() ; // draw bounding box
     
     } 
	board.saveEPS("fig2.eps");
}

int main( int argc, char** argv )
{
	srand(time(NULL));
	if (argc < 2) {
		trace.error() << "Please provide filename" << endl;
		return EXIT_FAILURE;
	}
	const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
	string filename = argv[1];
	float increment = atof(argv[2]);
	float noise = 0.0;
	if (argc >= 4)
		noise = atof(argv[3]);
	typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3D;
	typedef ImageFromSet<Image3D> SetConverter;
	typedef Pencil<Z3i::Point, MSTTangent<Z3i::Point>, PointVector<3, double> > Pencil; 
	int range = 200;
	int pitch =  20;
	int radius = 10;
	
	vector<PointVector<3, double>> curve;
//	createContinuousLogarithmicCurve(curve, 50, increment);
	createStraightLine(curve, 50, increment);
//	construct26ConnectedCurve(curve);
	set<PointVector<3,double>> vectorPoints;
	createVolumeFromCurve(curve, vectorPoints, 5);
//	thinVolume<Pencil>(curve, vectorPoints, 20.0);
//	drawDeformedCylinder(vectorPoints, 50, 5, increment);
//  createStraightLine(curve, 50);
	Z3i::Domain domain(Z3i::Point(-100,-100,-100), Z3i::Point(100, 300, 300));


	//createHelixCurve(vectorPoints, range, radius, pitch, increment);
//	drawCircle(vectorPoints, 50.0, 0., 0., 0., increment);
//	createSyntheticAirwayTree(vectorPoints, 5, 100, 0, 0, {0,0,0}, increment);
	Image3D anImage3D(domain);

//	create2DNaive();

/*	DigitalSet set(domain);
	for (auto it = vectorPoints.begin(), itE = vectorPoints.end(); it != itE; ++it) {
		set.insert(*it);
	}
	DigitalSet set2 = addNoise(set, noise);*/
	imageFromRangeAndValue(vectorPoints.begin(), vectorPoints.end(), anImage3D, 150);
    VolWriter<Image3D>::exportVol(examplesPath + filename, anImage3D);
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );


	trace.info() << "saved" << endl;
}
///////////////////////////////////////////////////////////////////////////////
