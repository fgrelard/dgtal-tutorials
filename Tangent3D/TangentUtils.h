#ifndef TANGENT_UTILS_H
#define TANGENT_UTILS_H

#include <vector>

namespace TangentUtils {
	double triangle(double x);
	double lambda(double x);

	template <typename Pencil, typename Iterator, typename Segmentation>
	std::vector<Pencil> computeTangentialCover(Iterator itB, Iterator itE,
											   const Segmentation& s);

}





/* Definition of the computation of tangential cover */
template <typename Pencil, typename Iterator, typename Segmentation>
std::vector<Pencil> TangentUtils::computeTangentialCover(Iterator itB, Iterator itE,
										   const Segmentation& s) {
	typedef typename Segmentation::SegmentComputerIterator SegmentComputerIterator;
	typedef typename Segmentation::SegmentComputer DSS;
	typedef typename DSS::Vector3d Point3D;
	typedef typename Pencil::T Tangent;
	typedef typename Pencil::Vector3d Vector3D;
	std::vector<Pencil> pencils;
	for (; itB != itE; ++itB) {
		//Pencil of tangents initialized, but used further
		Pencil pencil(*itB);
		std::vector<Tangent> tangents;
		for (SegmentComputerIterator sitB = s.begin(), sitE = s.end(); sitB != sitE; ++sitB) {
			//Size of the tangent in number of indices
			int size = 0;
			//Position in the tangent (index). Helps defining the eccentricity
			int position = 0;
			Tangent tangent;
			bool found = false;
			for (auto sitPointB = sitB->begin(), sitPointE = sitB->end(); sitPointB != sitPointE; ++sitPointB) {
				//If the point we re looking at is the same as the one in one of the segments
				if (*itB == *sitPointB) {
				
					DSS currentSegment(*sitB);
					//Direction vector
					// Getting the vector defining the segment
					Vector3D v = *(currentSegment.end()-1) - *currentSegment.begin();
					tangent = Tangent(position, v);
					found = true;
				}
				size++;
				position++;
			}
			if (found) {
				tangent.computeEccentricity(size);
				tangents.push_back(tangent);
			}
		}
		pencil.compute_tangent(tangents);
		pencils.push_back(pencil);
	}
	return pencils;
}

#endif
