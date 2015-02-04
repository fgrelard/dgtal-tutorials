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
#include "DGtal/io/writers/VolWriter.h"
///////////////////////////////////////////////////////////////////////////////


using namespace std;
using namespace DGtal;

double INCREMENT = 0.01;
/**
 * Board: where to draw2
 * range: range in z axis to draw the curve
 **/
template <typename Point>
class PointDeletable : public Point {
	typedef typename Point::Component Scalar;
public:
	PointDeletable() : Point() {}
	PointDeletable(const Scalar& x, const Scalar& y, const Scalar& z) : Point(x,y,z) {}
	bool myMarkedToDelete = false;
};

typedef DigitalSetSelector<Z3i::Domain,  HIGH_ITER_DS + MEDIUM_DS>::Type DigitalSet;
typedef Z3i::Space Space;
typedef Z3i::KSpace KSpace;
typedef Z3i::Curve::ArrowsRange Range;
typedef PointDeletable<Z3i::Point> Point;
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

void add(const Point& p, vector<Point> & v) {
	if (find(v.begin(), v.end(), p) != v.end()) {
		return;
	} else {
		v.push_back(p);
	}
}

/**
 * This function finds if the 6 neighbour of a given point p lie on the same orientation
 * If there are more than 2 direct neighbours, then we don't have the same orientation
 * that is to say we dont have a curve
 * This function allows to determine points which can be removed
 **/
 
bool sameDirection6Neighbour(const vector<Point> & v, Point & p) {
	typedef Z3i::Space Space;
	typedef MetricAdjacency<Space, 1> MetricAdjacency;

	vector<Point> neighbours6;
	for (auto it = v.begin(), itE = v.end(); it != itE; ++it) {
		if (MetricAdjacency::isAdjacentTo(*it, p) && *it != p) {
		    neighbours6.push_back(*it);
		}
	}

	bool sameDirection = true;
	for (unsigned int i = 0; i < neighbours6.size() - 1; i++) {
		Point current = neighbours6[i];
		Point next = neighbours6[i+1];
		if (!current.myMarkedToDelete && !next.myMarkedToDelete && current[0] != next[0] && current[1] != next[1]) {
			sameDirection = false;
			p.myMarkedToDelete = true;
		}
	}
	return sameDirection;
}

/**
 * Allows to obtain a 26 connected curve
 **/
void construct26ConnectedCurve(vector<Point> & vectorPoints) {
	vectorPoints.erase(
		std::remove_if(
			vectorPoints.begin(), vectorPoints.end(), [&vectorPoints](Point& other) {
				//the points removed correspond to those having neighbours with different direction
				return (!sameDirection6Neighbour(vectorPoints, other));
				
			}
			),
		vectorPoints.end()
		);
}

void drawCircle(vector<Point> & vectorPoints, float radius, float cx, float cy, float z) {
    float x, y;
	double angle = 0.0;
    while (angle <= 2 * M_PI){
		x = radius * cos( angle );
		y = radius * sin( angle );
		add(Point((int)x+cx, (int)y+cy, (int)z), vectorPoints);
		angle += INCREMENT;
	}
	construct26ConnectedCurve(vectorPoints);

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
		Point p (i, i, 0);
	    add(p, points);
	}
}


double lengthHelixCurve(int range, int radius, int pitch) {
	return range * sqrt(pow(radius, 2) + pow(pitch, 2));
}


int main( int argc, char** argv )
{
	if (argc < 2) {
		trace.error() << "Please provide filename" << endl;
		return EXIT_FAILURE;
	}
	const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
	string filename = argv[1];
	INCREMENT = atof(argv[2]);
	typedef ImageContainerBySTLVector<Z3i::Domain, unsigned char> Image3D;
	typedef ImageFromSet<Image3D> SetConverter;

	int range = 200;
	int pitch =  20;
	int radius = 10;
	
	vector<Point> vectorPoints;
	
	Z3i::Domain domain(Z3i::Point(-100,-100,-100), Z3i::Point(100, 100, 100));
	//Z3i::Domain domain(Z3i::Point(-100,-100,-100), Z3i::Point(100,100,100));

	//createStraightLine(vectorPoints, 50);
	//createHelixCurve(vectorPoints, range, radius, pitch);
	drawCircle(vectorPoints, 50.0, 0., 0., 0.);
	Image3D anImage3D(domain);
	imageFromRangeAndValue(vectorPoints.begin(), vectorPoints.end(), anImage3D, 150);
    //GenericWriter<Image3D>::exportFile(examplesPath + filename, anImage3D);
    VolWriter<Image3D>::exportVol(examplesPath + filename, anImage3D);
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	/*viewer.show();
	for (vector<Point>::iterator it = vectorPoints.begin() ; it != vectorPoints.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
		}*/


	trace.info() << "saved" << endl;
	//Display3dfactory<Space, KSpace>::draw(board, anImage3D);
    //VolWriter<ImageThree>::exportVol("aFilename.vol", anImage3D);
    //anImage3D >> "test.vol"; anImage3D << Z3i::Point(1,1,1);
}
///////////////////////////////////////////////////////////////////////////////
