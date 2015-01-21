/** File allowing to create a curve in DGtal according to an equation **/
#include <iostream>
#include <fstream>
#include "DGtal/base/Common.h"
#include "DGtal/io/Display3DFactory.h"
// Shape construction
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/curves/estimation/RosenProffittLocalLengthEstimator.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
// Drawing
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"
///////////////////////////////////////////////////////////////////////////////

#define INCREMENT 0.01
using namespace std;
using namespace DGtal;

/**
 * Board: where to draw2
 * range: range in z axis to draw the curve
 **/


typedef DigitalSetSelector<Z3i::Domain,  HIGH_ITER_DS + MEDIUM_DS>::Type DigitalSet;
typedef Z3i::Space Space;
typedef Z3i::KSpace KSpace;
typedef Z3i::Curve::ArrowsRange Range;
typedef Z3i::Point Point;
typedef Range::ConstIterator ConstIterator;
typedef RosenProffittLocalLengthEstimator<Range::ConstIterator> LengthEstimator; 
typedef vector<Point>::const_iterator PointIterator;
typedef  DigitalSetBySTLSet<DGtal::HyperRectDomain<DGtal::SpaceND<3u, int> >, std::less<DGtal::PointVector<3u, int, boost::array<int, 3ul> > > >::Iterator Iterator;

void saveToPL(vector<Point> & vPoints, std::string filename) {
	ofstream s(filename.c_str());
	for (vector<Point>::iterator it = vPoints.begin(); it != vPoints.end(); ++it) {
		s << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << endl;
	}
	s.close();
}

void add(const Point & p, vector<Point> & v) {
	if (find(v.begin(), v.end(), p) != v.end()) {
		return;
	} else {
		v.push_back(p);
	}
}

void drawCircle(vector<Point> & vectorPoints, float radius, float cx, float cy, float z) {
    float x, y;
	double angle = 0.0;
    while (angle <= 2 * M_PI){
		x = radius * cos( angle );
		y = radius * sin( angle );
		add(Point(x+cx, y+cy, z), vectorPoints);
		angle += 0.01f;
	}
}

/*
 * Radius Winding : radius of the loop
 * RadiusSpiral: radius of the volume
 */
void createHelixCurve(vector<Point> & vectorPoints, int range, int radiusWinding, int radiusSpiral, int pitch) {
	for (float i=0.0f; i < range; i+=INCREMENT) {
		float centerx = radiusWinding * cos(i/pitch);
		float centery = radiusWinding * sin(i/pitch);
		drawCircle(vectorPoints, radiusSpiral, centerx, centery, i);
	}
}


/*
 * This function does not allow to set the spiral radius
 */
void createHelixCurve( vector<Point> &point, int range, int radius, int pitch) {
	for (float i=0.0; i < range; i+=INCREMENT) {
	    int x = radius * cos(i);
		int y = radius * sin(i);
		int z = pitch * i;
		Point p(x,y,z);
		add(p, point);
	}
}

void createStraightLine(vector<Point> & points, int range) {
	for (int i = 0; i < range; i++) {
		Point p(i, i, 0);
	    add(p, points);
	}
}


double lengthHelixCurve(int range, int radius, int pitch) {
	return range * sqrt(pow(radius, 2) + pow(pitch, 2));
}

double estimateLength(vector<Point> v) {
    GridCurve<Z3i::KSpace> c;
	c.initFromVector(v);
	Range range = c.getArrowsRange();
	ConstIterator it = range.begin();
	ConstIterator itE = range.end();

//trace.info() << (*it).first << endl;
	

	LengthEstimator length;
	length.init(1, it, itE, true);
    return length.eval();
}

//template <typename KSPace, typename Space>
double estimateLength(const vector<Point> & points, /*Viewer3D<Space, KSPace> & viewer,*/ bool display=false) {
	typedef typename PointIterator::value_type Point3d;
	typedef StandardDSS6Computer<PointIterator,int,18> SegmentComputer;
	typedef SaturatedSegmentation<SegmentComputer> Decomposition;
	typedef typename Decomposition::SegmentComputerIterator SegmentComputerIterator;
	HueShadeColorMap<int> cmap_hue( 0, 3, 1 );
	SegmentComputer algo;
	Decomposition decomp(points.begin(), points.end(), algo);
	double totalEstimatedLength= 0.0;
//	viewer << SetMode3D( algo.className(), "BoundingBox" );
	int c = 0;
	for (SegmentComputerIterator i = decomp.begin(); i != decomp.end(); ++i) {
		SegmentComputerIterator tmp = i;
		SegmentComputer ms3d(*i);
		Point3d first = *ms3d.begin();
		Point3d last = *(ms3d.end() -1);
/*		trace.info() << "-  MS3D,"
					 << " [" << first[ 0 ] << "," << first[ 1 ] << ","<< first[ 2 ] << "]"
					 << "->[" << last[ 0 ] << "," << last[ 1 ] << ","<< last[ 2 ] << "]" << endl;*/
		double distanceDSS = sqrt(pow((last[0] - first[0]), 2) + pow((last[1] - first[1]), 2) + pow((last[2] - first[2]), 2));
		totalEstimatedLength += distanceDSS;
		Color color3d = cmap_hue(c);
		if (display) {
			//viewer << CustomColors3D(color3d, color3d) << ms3d;
		}
		c++;
	}
	return totalEstimatedLength;
}

int main( int argc, char** argv )
{
	if (argc < 2) {
		trace.error() << "Please provide filename" << endl;
		return EXIT_FAILURE;
	}
	const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
	string filename = argv[1];
	typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3D;
	typedef ImageFromSet<Image3D> SetConverter;

	int range = 200;
	int pitch =  20;
	int radius = 10;
	
	vector<Point> vectorPoints;
	
	Z3i::Domain domain(Z3i::Point(-radius,-radius,-pitch*range), Z3i::Point(radius,radius,pitch*range));
	//Z3i::Domain domain(Z3i::Point(-100,-100,-100), Z3i::Point(100,100,100));

	//createStraightLine(vectorPoints, 50);
    createHelixCurve(vectorPoints, range, radius, pitch);
	Image3D anImage3D(domain);
	imageFromRangeAndValue(vectorPoints.begin(), vectorPoints.end(), anImage3D, 150);
    //GenericWriter<Image3D>::exportFile(examplesPath + filename, anImage3D);
	saveToPL(vectorPoints, examplesPath + filename);
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	/*viewer.show();
	for (vector<Point>::iterator it = vectorPoints.begin() ; it != vectorPoints.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
		}*/


	trace.info() << "saved" << endl;
	trace.info() << "Expected= " << lengthHelixCurve(range, radius, pitch) << endl;
	trace.info() << "Measured= " << estimateLength(vectorPoints, /*viewer,*/ true) << endl;
	//Display3dfactory<Space, KSpace>::draw(board, anImage3D);
    //VolWriter<ImageThree>::exportVol("aFilename.vol", anImage3D);
    //anImage3D >> "test.vol"; anImage3D << Z3i::Point(1,1,1);
}
///////////////////////////////////////////////////////////////////////////////
