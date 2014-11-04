#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/boards/Board2D.h"


#include "DGtal/io/readers/PGMReader.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/Display2DFactory.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/geometry/volumes/distance/DistanceTransformation.h"
#include "DGtal/images/IntervalForegroundPredicate.h"

#include "DGtal/topology/helpers/Surfaces.h"
#include "DGtal/geometry/curves/estimation/DSSLengthEstimator.h"
#include <iostream>


const std::string examplesPath = "/home/florent/bin/DGtal/examples/samples/";

void tutoGridCurve() {
  DGtal::Z2i::Curve c;

  std::string square = examplesPath + "smallSquare.dat";

  std::fstream inputStream;

  //read file
  inputStream.open (square.c_str(), std::ios::in);

  //line drawn with coordinates x y
  c.initFromVectorStream(inputStream);
  inputStream.close(); 

  //Allows to display 2D objects
  DGtal::Board2D board; 
 
  board << c;
  DGtal::Z2i::Curve::OuterPointsRange r1 = c.getOuterPointsRange(); 
  board << r1; 
  board.saveEPS("gridCurve.eps");
}

void distanceTransformation() {
  
  /** Read a file **/
  typedef DGtal::ImageContainerBySTLVector< DGtal::Z2i::Domain, unsigned char> Image;
  typedef DGtal::GrayscaleColorMap<unsigned char> Gray;
  std::string filename =  examplesPath + "contourS.pgm";
  Image image = DGtal::PGMReader<Image>::importPGM(filename); 
  DGtal::trace.info() << "Imported image: "<<image<<std::endl;
  

  /** Saving domain and image **/
  DGtal::Board2D aBoard;
  aBoard << image.domain();  
  aBoard.saveSVG("imageDomainTuto.svg");
  aBoard.clear();
  DGtal::Display2DFactory::drawImage<Gray>(aBoard, image, (unsigned char)0, (unsigned char)255);
  aBoard.saveEPS("imageDomainTuto2.eps");

  /** Creating binarization and euclidean DT **/
  typedef DGtal::functors::IntervalForegroundPredicate<Image> Binarizer;
  //Threshold to 135
  Binarizer b(image,1, 135); 
  typedef DGtal::DistanceTransformation<DGtal::Z2i::Space, Binarizer, DGtal::Z2i::L2Metric> DTL2;
  DTL2 dt(&image.domain(),&b, &DGtal::Z2i::l2Metric );

  DTL2::Value maxDT = (*std::max_element(dt.constRange().begin(), 
                                         dt.constRange().end()));
  typedef DGtal::HueShadeColorMap<DTL2::Value,2> HueTwice;
  aBoard.clear();
  DGtal::Display2DFactory::drawImage<HueTwice>(aBoard, dt, (DTL2::Value)0, 
					       (DTL2::Value)maxDT);
  aBoard.saveEPS("imageDomainTuto3.eps");

}

void dssLength() {
  typedef DGtal::ImageContainerBySTLVector< DGtal::Z2i::Domain, unsigned char> Image;
  std::string filename = examplesPath + "contourS.pgm";
  Image image = DGtal::PGMReader<Image>::importPGM(filename);

  typedef DGtal::functors::IntervalThresholder<Image::Value> Binarizer;
  Binarizer b(1, 135);
  DGtal::functors::PointFunctorPredicate<Image, Binarizer> predicate(image, b);

  DGtal::Z2i::KSpace ks;
  ks.init( image.domain().lowerBound(), image.domain().upperBound(), true );
  DGtal::SurfelAdjacency<2> sAdj(true);

  std::vector< std::vector< DGtal::Z2i::SCell > > contours;
  DGtal::Surfaces<DGtal::Z2i::KSpace>::extractAll2DSCellContours(contours, ks, sAdj, predicate);

  DGtal::Z2i::Curve c;
  c.initFromSCellsVector( contours.at(1) );
  
  typedef DGtal::Z2i::Curve::PointsRange Range;
  Range r = c.getPointsRange();

  DGtal::DSSLengthEstimator< Range::ConstCirculator > DSSlength;
  DSSlength.init(1, r.c(), r.c());
  DGtal::trace.info() << "Length: " << DSSlength.eval() << std::endl;
  
}


void viewer3D() {
}

int main(int argc, char** argv) {

  //tutoGridCurve();
  //distanceTransformation();
  //dssLength();
  viewer3D();
  return 0;
}
