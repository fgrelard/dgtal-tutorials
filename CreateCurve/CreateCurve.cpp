/** File allowing to create a curve in DGtal according to an equation **/
#include <iostream>
#include "DGtal/base/Common.h"

// Shape construction
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/io/boards/Board3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/implicit/ImplicitBall.h"
#include "DGtal/shapes/EuclideanShapesDecorator.h"

// Drawing
#include "DGtal/io/boards/Board3D.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

/**
 * Board: where to draw
 * range: range in z axis to draw the curve
 **/

void createSpiralCurve(Board3D<> &board, int range) {
	for (int i; i < range; i++) {
		int x = cos(i);
		int y = sin(i);
		int z = i;
		board << Z3i::Point(x, y, z);
	}
}

int main( int argc, char** argv )
{
	if (argc < 2) {
		trace.error() << "Please provide filename" << endl;
		return EXIT_FAILURE;
	}
	const string examplesPath = "/home/florent/bin/DGtal/examples/samples/";
	string filename = argv[2];
	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	Board3D<> board;
	createSpiralCurve(board, 100);
	board.saveOBJ(examplesPath + filename, false);
}
///////////////////////////////////////////////////////////////////////////////
