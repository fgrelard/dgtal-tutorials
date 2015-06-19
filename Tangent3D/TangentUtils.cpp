#include <math.h>
#include "geometry/TangentUtils.h"

double TangentUtils::triangle(double x) {
	if (x >= 0 && x <= 0.5) {
		return (2 * x);
	} else if (x > 0.5 && x <= 1) {
		return ((1 - x) * 2);
	} else {
		return 0.;
	}
}

double TangentUtils::lambda(double x) {
	return 64*(pow(-x, 6) + pow(3 * x, 5) + pow(-3 * x, 4) + pow(x, 3));
}


