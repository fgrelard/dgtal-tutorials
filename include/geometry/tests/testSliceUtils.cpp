#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageSelector.h"
#include "../SliceUtils.h"

using namespace DGtal;
using namespace std;

void testSliceFromPlane() {
	typedef ImageSelector<Z3i::Domain, unsigned char>::Type Image;
	typedef DGtal::ConstImageAdapter<Image,Z2i::Domain,DGtal::functors::Point2DEmbedderIn3D<DGtal::Z3i::Domain>, Image::Value, DGtal::functors::Identity> ImageAdapterExtractor;

	Image image = GenericReader<Image>::import("/home/florent/test_img/cylinder.vol");

	Z3i::Point origin(0,0,1);
	Z3i::RealPoint normal(0,0,1);

	ImageAdapterExtractor extractedImage = SliceUtils::sliceFromPlane<ImageAdapterExtractor>(normal, origin, image, 100);
	GenericWriter<ImageAdapterExtractor>::exportFile("slice.pgm", extractedImage);	 
}


int main() {
//	testSliceFromPlane();
	return 0;
}
