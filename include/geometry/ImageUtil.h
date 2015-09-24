#ifndef IMAGE_UTIL_H
#define IMAGE_UTIL_H

namespace ImageUtil {
	template <typename ImageOut, typename ImageIn>
	ImageOut convertImage(const ImageIn& imageIn, int factor = 1);
}

template <typename ImageOut, typename ImageIn>
ImageOut ImageUtil::convertImage(const ImageIn& imageIn, int factor) {
	ImageOut imageOut(imageIn.domain());

	for (auto it = imageIn.domain().begin(), ite = imageIn.domain().end();
		 it != ite; ++it) {
		imageOut.setValue(*it, imageIn(*it) * factor);
	}
	return imageOut;
}

#endif
