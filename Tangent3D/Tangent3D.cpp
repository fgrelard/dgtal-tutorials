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
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"

// Local includes
#include "Pencil.h"
#include "MSTTangent.h"
#include "TangentUtils.h"
#include "SliceUtils.h"

using namespace DGtal;
using namespace std;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GreedySegmentation
///////////////////////////////////////////////////////////////////////////////


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

template <typename Pencil>
void visualize(const vector<Pencil> & tangents, Viewer3D<> & viewer) {
	typedef DGtal::PointVector<3, double> Vector3D;
	typedef DGtal::PointVector<3, int> Point3D;
	int size_factor = 0;
	int plane_factor = 3;
	for (auto it = tangents.begin(), itE = tangents.end(); it != itE; ++it) {
		if (!it->isUndefined()) {
			// SegmentComputer currentSegmentComputer(*i); 
			// viewer << SetMode3D(currentSegmentComputer.className(), "BoundingBox");
			viewer << CustomColors3D(Color::Red, Color::Red);
			//viewer << currentSegmentComputer;
			Point3D point = it->getPoint();
			Vector3D tangent = it->getTangent();

			viewer.addLine(point - (tangent * size_factor), point + (tangent * size_factor), 0.1);
			vector<Vector3D> plane = computePlaneFromNormalVector<Vector3D>(it->getTangent());
			Color planeColor(0, 0, 255, 128);
			//viewer << CustomColors3D(planeColor, planeColor);

			
			Vector3D p1 = (Vector3D) point + plane[0] * plane_factor;
			Vector3D p2 = (Vector3D) point + plane[1] * plane_factor;
			Vector3D p3 = (Vector3D) point + plane[2] * plane_factor;
			Vector3D p4 = (Vector3D) point + plane[3] * plane_factor;
			viewer.addQuad(p1, p2, p3, p4);
		}
	}
}


/**
 * saturated segmentation of a (sub)range
 */
template <typename Pencil, typename Iterator>
vector<Pencil> orthogonalPlanesWithTangents(Iterator itb, Iterator ite, Viewer3D<> & viewer)
{
	typedef StandardDSS6Computer<Iterator,int,8> SegmentComputer;  
	typedef SaturatedSegmentation<SegmentComputer> Segmentation;

	SegmentComputer algo;
	Segmentation s(itb, ite, algo);
	s.setMode("MostCentered++");
	typename Segmentation::SegmentComputerIterator i = s.begin();
	typename Segmentation::SegmentComputerIterator end = s.end();
	
	
	for (; i != end; ++i) {
		SegmentComputer currentSegmentComputer(*i); 
		viewer << SetMode3D(currentSegmentComputer.className(), "BoundingBox");
		viewer << CustomColors3D(Color::Green, Color::Green);
	   	viewer << currentSegmentComputer;
	}

	
	vector<Pencil> tangents = TangentUtils::computeTangentialCover<Pencil>(itb, ite, s);
	//visualize(tangents, viewer);
    return tangents;
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
		("input2,v", po::value<std::string>(), "vol file (corresponding volume)")
		("output,o", po::value<std::string>(), "sliced vol file with orthogonal planes")
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
	typedef DGtal::PointVector<3, double> Vector3D;
	typedef DGtal::PointVector<3, int> Point3D;
	typedef MSTTangent<Point3D> Tangent;
	typedef Pencil<Vector3D, Tangent, Vector3D> Pencil;
	
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
	vector<Pencil> tangents = orthogonalPlanesWithTangents<Pencil>(vPoints.begin(), vPoints.end(), viewer);
//	visualize(tangents, viewer);
	if (vm.count("output") && vm.count("input2")) {
		string inputFileName2 = vm["input2"].as<std::string>();
		Image volume = GenericReader<Image>::import(inputFileName2);
		string outName = vm["output"].as<std::string>();
		SliceUtils::slicesFromPlanes(viewer, tangents, volume, outName);
	}
	
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	for (auto it = vPoints.begin(); it != vPoints.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
	}
	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();

	return 0;
}
