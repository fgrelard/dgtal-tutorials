

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

/**
 * @file meshFromOFF.cpp
 * @ingroup Examples
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2012/07/05
 *
 * An example file named meshFromOFF.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////


//! [includeImportOFF]
//!


#include "DGtal/io/readers/MeshReader.h"

#include <QtGui/qapplication.h>
#include "DGtal/io/Display3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
//! [includeImportOFF]
#include "DGtal/base/Common.h"
#include "DGtal/io/Color.h"
///////////////////////////////////////////////////////////////////////////////
using namespace std;
using namespace DGtal;
///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
	const std::string examplesPath = "/home/florent/bin/DGtal/examples/samples/"
  QApplication application(argc,argv);
  Viewer3D<> viewer;
  viewer.show();     
  //! [ImportOFFfile]
  std::string inputFilename = examplesPath + "tref.off";   
  // Since the input points are not necessary integers we use the PointD3D from Display3D.
  Mesh<Viewer3D<>::RealPoint> anImportedMesh;
  anImportedMesh << inputFilename;
  //! [ImportOFFfile]
  trace.info()<< "importating done..."<< endl;
  //! [  viewer.setLineColor(DGtal::Color(150,0,0,254));
  viewer << anImportedMesh;
  viewer << Viewer3D<>::updateDisplay;
  //! [displayOFFfile]
  return application.exec();
}
///////////////////////////////////////////////////////////////////////////////