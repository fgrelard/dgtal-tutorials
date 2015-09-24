#ifndef __MORPHOMATHS__
#define __MORPHOMATHS__

#include <limits>
#include "DGtal/base/Common.h"
#include "DGtal/base/BasicFunctors.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/kernel/domains/DomainPredicate.h"

namespace Morphomaths {	
	template <typename Image, typename Domain>
	bool process(const Image& image, const Domain& domain, int x, int y, int nx, int ny);

	template <typename Image>
	Image constructOnePxBorderImage(const Image& image);

	template <typename Image>
	Image erosion(const Image& image, int size);

	template <typename Image>
	Image dilation(const Image& image, int size);
	
	template <typename Image>
	Image open(const Image & image, int size);

	template <typename Image>
	Image close(const Image& image, int size);

}


template <typename Image, typename Domain>
bool Morphomaths::process(const Image& image, const Domain& domain, int x, int y, int nx, int ny) {
	typedef typename Image::Point Point;
	
	int valueMin = std::numeric_limits<int>::max();
    int valueMax = 0;
	for (int i = -nx; i <= nx; i++) {
		for (int j = -ny; j <= ny; j++) {
			Point current(x+i,y+j);
			if (domain.isInside(current)) {
				int value = image(current);
				if ((i != 0 || j != 0)) {
					if (value < valueMin) valueMin = value;
					if (value > valueMax) valueMax = value;
				}
			}
		}
	}
	return (valueMin != valueMax);
}

template <typename Image>
Image Morphomaths::constructOnePxBorderImage(const Image& image) {
	typedef typename Image::Domain Domain;
	typedef typename Image::Point Point;
	
	Domain domain = image.domain();
	
	Image toReturn(Domain(domain.lowerBound() - Point::diagonal(), domain.upperBound() + Point::diagonal()));
	for (auto it = toReturn.domain().begin(), ite = toReturn.domain().end(); it != ite; ++it) {
		if (domain.isInside(*it))
			toReturn.setValue(*it, image(*it));
		else
			toReturn.setValue(*it, 0);
	}

	return toReturn;
}

template <typename Image>
Image Morphomaths::erosion(const Image& image, int size) {
	typedef typename Image::Domain Domain;
	typedef typename Image::Point Point;

	Domain domain = image.domain();
	Image toReturn = constructOnePxBorderImage(image);
	Image toWork = toReturn;
	Point upper = domain.upperBound(), lower = domain.lowerBound();

	int width = upper[0] - lower[0]+1;
	int height = upper[1] - lower[1]+1;

	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			bool shouldBeEroded = process(toWork, toWork.domain(), i, j, size, size);
			if (shouldBeEroded) toReturn.setValue(Point(i,j), 0);
		}
	}
	Image out(domain);
	for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
		out.setValue(*it, toReturn(*it));
	}
	return out;
}

template <typename Image>
Image Morphomaths::dilation(const Image& image, int size) {
	typedef typename Image::Domain Domain;
	typedef typename Image::Point Point;	
	typedef DGtal::functors::NotPointPredicate<Image> BackgroundPredicate; 
	
	
	Image toReturn = constructOnePxBorderImage(image);
	Image toWork = toReturn;
	BackgroundPredicate backgroundPredicate( toWork );
	
	Domain domain = image.domain();
	Point upper = domain.upperBound(), lower = domain.lowerBound();
	int width = upper[0] - lower[0] + 1;
	int height = upper[1] - lower[1] +  1;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			bool shouldBeDilated = process(backgroundPredicate, toWork.domain(), i, j, size, size);
			if (shouldBeDilated) toReturn.setValue(Point(i,j), 1);
		}
	}
	Image out(domain);
	for (auto it = domain.begin(), ite = domain.end(); it != ite; ++it) {
		out.setValue(*it, toReturn(*it));
	}
	return out;
}

template <typename Image>
Image Morphomaths::open(const Image & image, int size) {
	Image eros = erosion(image, size);
	Image dilat = dilation(eros, size);
	return dilat;
}

template <typename Image>
Image Morphomaths::close(const Image& image, int size) {
	Image dilat = dilation(image, size);
	Image eros = erosion(dilat, size);
	return dilat;
}

#endif
