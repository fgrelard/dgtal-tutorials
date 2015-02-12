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
#include <Eigen/Geometry>
#include <Eigen/Dense>

using namespace std;
using namespace DGtal;

double INCREMENT = 0.01;
template <typename Point>
double euclideanDistance(Point, Point);

template <typename Point>
class PointDeletable : public Point {
	typedef typename Point::Component Scalar;
public:
	PointDeletable() : Point() {}
	PointDeletable(const Scalar& x, const Scalar& y, const Scalar& z) : Point(x,y,z) {}
	bool myMarkedToDelete = false;
};

template <typename Point>
class Ball {
public:
	Ball() : myRadius(0), myCenter({0,0,0}) {}
	Ball(const Point& center, int radius) : myCenter(center), myRadius(radius) {}
	bool contains(const Point& point) const {return euclideanDistance(point, myCenter) <= myRadius;}
	vector<Point> pointsInBall() const {
		vector<Point> points;
		for (int i = -myRadius + myCenter[0], iend = myRadius + myCenter[0] + 1; i < iend; i++) {
			for (int j = -myRadius + myCenter[1], jend = myRadius + myCenter[1] + 1; j < jend; j++) {
				for (int k = -myRadius + myCenter[2], kend = myRadius + myCenter[2] + 1; k < kend; k++) {
					Point p(i, j, k);
					if (contains(p)) {
						points.push_back(p);
					}
				}
			}
		}
		return points;
	};
	bool operator!=(const Ball & other) const {return (myCenter != other.myCenter || myRadius != other.myRadius);}
private:
	Point myCenter;
	int myRadius;
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

template <typename Point>
double euclideanDistance(Point p, Point other) {
	return sqrt(pow( (p[0] - other[0]), 2) + pow( (p[1] - other[1]), 2) + pow( (p[2] - other[2]), 2) );
}

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

double euclideanDistance(float x1, float y1, float x2, float y2) {
	return sqrt( pow( (x2 - x1), 2) + pow( (y2 - y1), 2));
}

void drawDisk(Eigen::Matrix<float, Eigen::Dynamic, 4> & m, float radius, float cx, float cy, float cz, long int &row) {
	for (int x = cx - radius, xend = cx + radius + 1; x < xend; x++) {
		for (int y = cy - radius, yend = cy + radius +  1; y < yend; y++) {
			if (euclideanDistance(x, y, cx, cy) <= radius) {
				m(row, 0) = x;
				m(row, 1) = y;
				m(row, 2) = cz;
				m(row, 3) = 0;
				row++;
				m.conservativeResize(row+1, 4);
			}
		}
	}
}

void drawCylinder(Eigen::Matrix<float, Eigen::Dynamic, 4> & m, int length, int radius, float cz, float rotationAngle) {
	int progress = 0;
	long int row = 0;
    while (progress < length) {
		drawDisk(m, radius, 0, 0, cz+progress, row);
		progress++;
	}
	//Deletes
	m.conservativeResize(row, 4);
	//Rotation with rotationAngle
	Eigen::Affine3f rot(Eigen::AngleAxisf(rotationAngle, Eigen::Vector3f::UnitX()));
	if (rotationAngle != 0.0) {
		m = m * rot.matrix();
	}
}


/*
 * Radius  Winding : radius of the loop
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

void createSyntheticAirwayTree(vector<Point> & v, int branchNumber, int lengthTrachea, int z, float rotationAngle, Point firstPoint) {
	if (branchNumber == 0) return;
	int radius = lengthTrachea / 6;
	cout << radius << " " << rotationAngle << endl;
	Eigen::Matrix<float, Eigen::Dynamic, 4> *matrix = new Eigen::Matrix<float, Eigen::Dynamic, 4>(1, 4);
	//Creates a rotated cylinder
	drawCylinder(*matrix, lengthTrachea, radius, z, rotationAngle);
	//Filling the vector with points from matrix
	for (int i = 0; i < matrix->innerSize(); i++) {
		if ((*matrix)(i, 0) != 0 || (*matrix)(i, 1) != 0 || (*matrix)(i, 2) != 0) {
			int posx = (*matrix)(i, 0);
			//Translation after rotation
			int posy = (*matrix)(i, 1) - (z - ((radius * 2) / sqrt(2))) * sin(rotationAngle) + firstPoint[1];
			int posz = (*matrix)(i, 2) - (z - ((radius * 2) / sqrt(2))) * cos(rotationAngle) + firstPoint[2];
			add(Point(posx, posy, posz), v);
		}
	}
	z += lengthTrachea;
	//Determining initial starting point for new branches
	firstPoint += Point(0, lengthTrachea * sin(rotationAngle), lengthTrachea * cos(rotationAngle));
	//matrix.resize(0, 0);
	matrix->resize(0,0);
	delete matrix;
	createSyntheticAirwayTree(v, branchNumber - 1, lengthTrachea * 0.6, z, rotationAngle + 0.2 * M_PI, firstPoint);
	createSyntheticAirwayTree(v, branchNumber - 1, lengthTrachea * 0.6, z, rotationAngle - 0.2 * M_PI, firstPoint);
}

void createLogarithmicCurve(vector<Point> & curve, int range) {
	for (float i = 1; i < range; i+=INCREMENT) {
		float x = i;
		float y = i;
		float z = 20*log(i);
		add(Point((int)x, (int)y, (int)z), curve);
	}
}

void createVolumeFromCurve(const vector<Point> & curve, set<Point> & volume, int ballRadius) {
	vector<Ball<Point> > ballVector;
	for (Point point : curve) {
		ballVector.push_back(Ball<Point>(point, ballRadius));
	}
	for (const Ball<Point>& current : ballVector) {
		for (Point& point : current.pointsInBall()) {
			for (const Ball<Point>& other : ballVector) {
				if (other != current && !other.contains(point)) {
					volume.insert(point);
				}
			}
		}
	}
	
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
	
	//set<Point> vectorPoints;
	vector<Point> curve;
	//createLogarithmicCurve(curve, 50);
	//createVolumeFromCurve(curve, vectorPoints, 20);

	vector<Point> vectorPoints;
	Z3i::Domain domain(Z3i::Point(-100,-100,-100), Z3i::Point(100, 100, 300));


	//createHelixCurve(vectorPoints, range, radius, pitch);
//	drawCircle(vectorPoints, 50.0, 0., 0., 0.);
	createSyntheticAirwayTree(vectorPoints, 5, 100, 0, 0, {0,0,0});
	Image3D anImage3D(domain);
	imageFromRangeAndValue(vectorPoints.begin(), vectorPoints.end(), anImage3D, 150);
    VolWriter<Image3D>::exportVol(examplesPath + filename, anImage3D);
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );


	trace.info() << "saved" << endl;
}
///////////////////////////////////////////////////////////////////////////////
