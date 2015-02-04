/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
//LICENSE-END
/**
 * @file testSegmentation.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 *
 * @date 2011/07/22
 *
 * This file is part of the DGtal library
 */

/**
 * Description of testSegmentation <p>
 * Aim: simple test of \ref GreedySegmentation and \ref SaturatedSegmentation
 */

#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <iterator>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/base/Common.h"
#include "DGtal/base/Exceptions.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/base/Circulator.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/VolReader.h" 
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/shapes/parametric/Ball2D.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/shapes/Shapes.h"

using namespace DGtal;
using namespace std;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GreedySegmentation
///////////////////////////////////////////////////////////////////////////////

double triangle(double);
double lambda(double);

template <typename Vector3D>
class MSTTangent {
public:
	typedef Vector3D Vector;
	MSTTangent() : position(0), v(Vector()) {}
	MSTTangent(int _position, Vector _v) : position(_position), v(_v) {}
	void computeEccentricity(float size) { eccentricity = position / size; }
private:
	int position;
public:
	Vector v;
	float eccentricity;
};

template <typename Point, typename Tangent, typename Vector>
class Pencil {
public:
	typedef Tangent T;
	typedef Vector Vector3d;
public:
	Pencil(Point _point) : point(_point) {}
	Vector getTangent() { return tangent; }
	Point getPoint() { return point; }
	bool isUndefined() { return undefined; }
	void compute_tangent(vector<Tangent> tangents) {
		Vector v;
	    double weighted_eccentricity = 0.;
		
		for (auto it = tangents.begin(), itE = tangents.end(); it != itE; ++it) {
			float weight = triangle(it->eccentricity);
			Vector direction = (Vector) it->v / (float) it->v.norm();
		    direction *= weight; 
			v += direction;	
			weighted_eccentricity += weight;
		}
	
		if (weighted_eccentricity != 0.)
			tangent = v / weighted_eccentricity;
		else if (v == v) { // checking if nan
			tangent = v;
		}
		if (tangent == Vector()) 
			undefined = true;
		
	}
private:
    Vector tangent;
	Point point;
	bool undefined = false;
};


template <typename Pencil, typename Iterator, typename Segmentation>
vector<Pencil> computeTangentialCover(Iterator itB, Iterator itE,
							const Segmentation& s) {
	typedef typename Segmentation::SegmentComputerIterator SegmentComputerIterator;
	typedef typename Segmentation::SegmentComputer DSS;
	typedef typename DSS::Vector3d Point3D;
	typedef typename Pencil::T Tangent;
	typedef typename Pencil::Vector3d Vector3D;
	vector<Pencil> pencils;
	for (; itB != itE; ++itB) {
		//Pencil of tangents initialized, but used further
		Pencil pencil(*itB);
		vector<Tangent> tangents;
		for (SegmentComputerIterator sitB = s.begin(), sitE = s.end(); sitB != sitE; ++sitB) {
			//Size of the tangent in number of indices
			int size = 0;
			//Position in the tangent (index). Helps defining the eccentricity
			int position = 0;
			Tangent tangent;
			bool found = false;
			for (auto sitPointB = sitB->begin(), sitPointE = sitB->end(); sitPointB != sitPointE; ++sitPointB) {
				//If the point we re looking at is the same as the one in one of the segments
				if (*itB == *sitPointB) {
				
					DSS currentSegment(*sitB);
					//Direction vector
					// Getting the vector defining the segment
					Vector3D v = *(currentSegment.end()-1) - *currentSegment.begin();
					tangent = Tangent(position, v);
					found = true;
				}
				size++;
				position++;
			}
			if (found) {
				tangent.computeEccentricity(size);
				tangents.push_back(tangent);
			}
		}
		pencil.compute_tangent(tangents);
		pencils.push_back(pencil);
	}
	return pencils;
}

float deduceCoordinateFromEquationThreeUnknown(float a, float b, float c, float first_value, float second_value) {
	return (-(a * first_value + b * second_value) / c); 
}

template <typename Point, typename Vector>
vector<Point> computePlaneFromNormalVector(Vector normal) {
	double a = normal[0];
	double b = normal[1];
	double c = normal[2];
	vector<Point> fourPointsForPlane;

	float x, y, z;
	Point p1, p2, p3, p4;
	if (a != 0) {
		y = 1;
		z = 1;
		x = deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p1 = Vector(x, y, z).getNormalized();
		y = -1;
		x = deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p2 = Vector(x, y, z).getNormalized();
		z = -1;
		x = deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p3 = Vector(x, y, z).getNormalized();
		y = 1;
		x = deduceCoordinateFromEquationThreeUnknown(b, c, a, y, z);
		p4 = Vector(x, y, z).getNormalized();
	} else if (b != 0) {
		x = 1;
		z = 1;
		y = deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p1 = Vector(x, y, z).getNormalized();
		x = -1;
		y = deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p2 = Vector(x, y, z).getNormalized();
		z = -1;
		y = deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p3 = Vector(x, y, z).getNormalized();
		x = 1;
		y = deduceCoordinateFromEquationThreeUnknown(a, c, b, x, z);
		p4 = Vector(x, y, z).getNormalized();
	} else if (c != 0) {
		x = 1;
		y = 1;
		z = deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p1 = Vector(x, y, z).getNormalized();
		y = -1;
		z = deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p2 = Vector(x, y, z).getNormalized();
		x = -1;
		z = deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p3 = Vector(x, y, z).getNormalized();
		y = 1;
		z = deduceCoordinateFromEquationThreeUnknown(a, b, c, x, y);
		p4 = Vector(x, y, z).getNormalized();
	}
	fourPointsForPlane = {p1, p2, p3, p4};
    
	return fourPointsForPlane;
}


/**
 * saturated segmentation of a (sub)range
 */
template <typename Iterator>
void orthogonalPlanesWithTangents(Iterator itb, Iterator ite, Viewer3D<> & viewer)
{
	typedef StandardDSS6Computer<Iterator,int,8> SegmentComputer;  
	typedef SaturatedSegmentation<SegmentComputer> Segmentation;

	SegmentComputer algo;
	Segmentation s(itb, ite, algo);
	s.setMode("MostCentered++");
	typename Segmentation::SegmentComputerIterator i = s.begin();
	typename Segmentation::SegmentComputerIterator end = s.end();
	typedef DGtal::PointVector<3, double> Vector3D;
	typedef DGtal::PointVector<3, int> Point3D;
	typedef MSTTangent<Point3D> Tangent;
	typedef Pencil<Vector3D, Tangent, Vector3D> Pencil;
	
	/*for (; i != end; ++i) {
		SegmentComputer currentSegmentComputer(*i); 
		viewer << SetMode3D(currentSegmentComputer.className(), "BoundingBox");
		viewer << CustomColors3D(Color::Green, Color::Green);
		viewer << currentSegmentComputer;
	}*/

	
	vector<Pencil> tangents = computeTangentialCover<Pencil>(itb, ite, s);
	for (auto it = tangents.begin(), itE = tangents.end(); it != itE; ++it) {
		if (!it->isUndefined()) {
			// SegmentComputer currentSegmentComputer(*i); 
			// viewer << SetMode3D(currentSegmentComputer.className(), "BoundingBox");
			viewer << CustomColors3D(Color::Red, Color::Red);
			//viewer << currentSegmentComputer;
			Point3D point = it->getPoint();
			Vector3D tangent = it->getTangent();
			int size_factor = 0;
			viewer.addLine(point - (tangent * size_factor), point + (tangent * size_factor), 0.1);
			vector<Vector3D> plane = computePlaneFromNormalVector<Vector3D>(it->getTangent());
			Color planeColor(0, 0, 255, 128);
			viewer << CustomColors3D(planeColor, planeColor);

			int plane_factor = 3;
			Vector3D p1 = (Vector3D) point + plane[0] * plane_factor;
			Vector3D p2 = (Vector3D) point + plane[1] * plane_factor;
			Vector3D p3 = (Vector3D) point + plane[2] * plane_factor;
			Vector3D p4 = (Vector3D) point + plane[3] * plane_factor;
			viewer.addQuad(p1, p2, p3, p4);
		}
	}
}

double triangle(double x) {
	if (x >= 0 && x <= 0.5) {
		return (2 * x);
	} else if (x > 0.5 && x <= 1) {
		return ((1 - x) * 2);
	} else {
		return 0.;
	}
}

double lambda(double x) {
	return 64*(pow(-x, 6) + pow(3 * x, 5) + pow(-3 * x, 4) + pow(x, 3));
}

template <typename Point>
	vector<Point> createBall(int x, int y, float radius) {
	typedef typename Z3i::Space MySpace;
	typedef Ball2D<MySpace> Ball;
	typedef Z3i::Domain Domain;
    Ball ball(x, y, radius);
    GaussDigitizer<MySpace, Ball> digitizer;
	digitizer.attach(ball);
	digitizer.init( ball.getLowerBound() + Z3i::Vector(-1,-1,-1),
					ball.getUpperBound() + Z3i::Vector(1, 1, 1),
					1 );
	Z3i::DigitalSet dgtalSet(digitizer.getDomain());
	Shapes<Domain>::digitalShaper(dgtalSet, digitizer);
for (auto it = dgtalSet.begin(), itE = dgtalSet.end();
	 it != itE;
	 ++it) {
	cout << *it << endl;
}
	
}

/////////////////////////////////////////////////////////////////////////
//////////////// MAIN ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (skeleton)")
		("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
		("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
		; 

	bool parseOK=true;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, general_opt), vm);  
	} catch(const std::exception& ex){
		parseOK=false;
		trace.info()<< "Error checking program options: "<< ex.what()<< endl;
	}
	po::notify(vm);    
	if( !parseOK || vm.count("help")||argc<=1)
	{
		std::cout << "Usage: " << argv[0] << " [input]\n"
				  << "Display volume file as a voxel set by using QGLviewer"<< endl
				  << general_opt << "\n";
		return 0;
	}  
	if(!vm.count("input"))
	{
		trace.error() << " The file name was not defined" << endl;      
		return 0;
	}
	string inputFilename = vm["input"].as<std::string>();
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	
	typedef PointVector<3,int> Point;
	QApplication application(argc,argv);
	Viewer3D<> viewer;
	viewer.show();
	
	
	typedef Z3i::Space Space;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef Z3i::KSpace KSpace;
	KSpace ks;
	
	Image image = VolReader<Image>::importVol(inputFilename);

	//Extracts the first point belonging to the object
	Point p;	
	vector<Point> vPoints;
	for (auto it = image.domain().begin(), itE = image.domain().end(); it != itE; ++it) {
		if (image(*it) >= thresholdMin && image(*it) <= thresholdMax) {
		    p = *it;
		}
	}
	// We have to visit the direct neighbours in order to have a container with voxels
	// ordered sequentially by their connexity
	// Otherwise we have a point container with points which are not neighbours
	// and this impairs maximal segment recognition
	typedef MetricAdjacency<Space, 3> Graph;
	typedef DepthFirstVisitor<Graph, set<Point> > Visitor;
	typedef typename Visitor::Node MyNode;
	typedef GraphVisitorRange<Visitor> VisitorRange;
	Graph graph;
	Visitor visitor( graph, p );
	MyNode node;

    
	while ( !visitor.finished() ) 
	{
  		node = visitor.current();
		if ( image.domain().isInside(node.first) &&
			 image(node.first) >= thresholdMin &&
			 image(node.first) <= thresholdMax ) { //is inside domain
			vPoints.push_back(node.first);
			visitor.expand();
		}
		else
			visitor.ignore();
	}
	orthogonalPlanesWithTangents(vPoints.begin(), vPoints.end(), viewer);
	
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	for (auto it = vPoints.begin(); it != vPoints.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
	}
	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();

	return 0;
}
