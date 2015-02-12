#ifndef SLICE_UTILS_H
#define SLICE_UTILS_H

#include <iostream>
#include <vector>
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/images/ImageHelper.h"
//! [ExampleViewer3D2DImagesExtractImagesNonSliceHeader]
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/io/viewers/Viewer3D.h"

using namespace DGtal;

namespace SliceUtils {
	template <typename Pencil, typename Image>
	void slicesFromPlanes(Viewer3D<>&, const std::vector<Pencil> &, const Image&, std::string);
}

template <typename Pencil, typename Image>
void SliceUtils::slicesFromPlanes(Viewer3D<>& viewer, const std::vector<Pencil> & vectorPlanes, const Image& volume, std::string outFileName) {
	typedef Image Image3D;
	typedef typename Image3D::Value Value;
	//! [ExampleViewer3D2DImagesExtractImagesNonSliceType]
	typedef DGtal::ConstImageAdapter<Image3D,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>,Value, DGtal::functors::Identity> ImageAdapterExtractor;


	DGtal::functors::Identity idV;
  

  const int IMAGE_PATCH_WIDTH = 20;  
  // Setting the image domain of the resulting image to be displayed in 3D:
  DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
                                    DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH)); 
  //! [ExampleViewer3D2DImagesExtractImagesNonSliceParam]
  
  
  
  unsigned int sliceNumber = 0;
	
  for (auto it = vectorPlanes.begin(), itE = vectorPlanes.end(); it != itE; ++it) {
typename Pencil::Vector3d planeNormal = it->getTangent();

	  typename Pencil::P origin = it->getPoint();
	  DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain >  embedder(volume.domain(), origin, planeNormal, IMAGE_PATCH_WIDTH);
	  ImageAdapterExtractor extractedImage(volume, domainImage2D, embedder, idV);
	  std::string outName;
	  outName += outFileName + "_" + std::to_string(sliceNumber) + ".pgm";
	  //GenericWriter<ImageAdapterExtractor>::exportFile(outName, extractedImage);
	  sliceNumber++;
	  DGtal::Z2i::Domain domainImage2D (DGtal::Z2i::Point(0,0), 
                                    DGtal::Z2i::Point(IMAGE_PATCH_WIDTH, IMAGE_PATCH_WIDTH)); 
	  viewer << extractedImage;
	  viewer << DGtal::UpdateImage3DEmbedding<Z3i::Space, Z3i::KSpace>(sliceNumber, 
                                                                     embedder(Z2i::RealPoint(0,0)),
                                                                     embedder(Z2i::RealPoint(IMAGE_PATCH_WIDTH,0)),
                                                                     embedder(domainImage2D.upperBound()),
                                                                     embedder(Z2i::RealPoint(0, IMAGE_PATCH_WIDTH)));
	  if (sliceNumber >= 10) break;
  }
}

#endif
