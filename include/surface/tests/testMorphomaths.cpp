#include "../Morphomaths.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

using namespace DGtal;

void testErosion() {
	typedef ImageContainerBySTLVector<Z2i::Domain, bool> Image;
	typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> ImageOut;
	Image image = GenericReader<Image>::import("/home/florent/trash/slice_1.pgm");
	Image ero = Morphomaths::erosion(image, 1);

	ImageOut out(ero.domain());
	for (auto it = out.domain().begin(), ite = out.domain().end(); it != ite; ++it) {
		if (ero(*it) == 0) {
			out.setValue(*it, 0);
		}
		else {
			out.setValue(*it, 255);
		}
	}
	
	GenericWriter<ImageOut>::exportFile("/home/florent/trash/bero.pgm", out);
	trace.info() << "done " << std::endl;
}

 void testDilation() {
	typedef ImageContainerBySTLVector<Z2i::Domain, bool> Image;
	typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> ImageOut;
	Image image = GenericReader<Image>::import("/home/florent/trash/slice_1.pgm");
	Image ero = Morphomaths::dilation(image, 1);

	ImageOut out(ero.domain());
	for (auto it = out.domain().begin(), ite = out.domain().end(); it != ite; ++it) {
		if (ero(*it) == 0) {
			out.setValue(*it, 0);
		}
		else {
			out.setValue(*it, 255);
		}
	}
	
	GenericWriter<ImageOut>::exportFile("/home/florent/trash/bdil.pgm", out);
 }

void testOpen() {
	typedef ImageContainerBySTLVector<Z2i::Domain, bool> Image;
	typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> ImageOut;
	Image image = GenericReader<Image>::import("/home/florent/trash/slice_1.pgm");
	Image ero = Morphomaths::open(image, 1);

	ImageOut out(ero.domain());
	for (auto it = out.domain().begin(), ite = out.domain().end(); it != ite; ++it) {
		if (ero(*it) == 0) {
			out.setValue(*it, 0);
		}
		else {
			out.setValue(*it, 255);
		}
	}
	
	GenericWriter<ImageOut>::exportFile("/home/florent/trash/bopen.pgm", out);
}
 
	
void testClose() {
	 typedef ImageContainerBySTLVector<Z2i::Domain, bool> Image;
	typedef ImageContainerBySTLVector<Z2i::Domain, unsigned char> ImageOut;
	Image image = GenericReader<Image>::import("/home/florent/trash/slice_1.pgm");
	Image ero = Morphomaths::close(image, 1);

	ImageOut out(ero.domain());
	for (auto it = out.domain().begin(), ite = out.domain().end(); it != ite; ++it) {
		if (ero(*it) == 0) {
			out.setValue(*it, 0);
		}
		else {
			out.setValue(*it, 255);
		}
	}
	
	GenericWriter<ImageOut>::exportFile("/home/florent/trash/bclose.pgm", out);
 }

int main() {
	testErosion();
	testDilation();
	testClose();
	testOpen();
	return 0;
}
