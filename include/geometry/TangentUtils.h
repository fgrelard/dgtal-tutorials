#ifndef TANGENT_UTILS_H
#define TANGENT_UTILS_H

#include <vector>
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/base/Common.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/helpers/StdDefs.h"

using namespace DGtal;

namespace TangentUtils {
	double triangle(double x);
	double lambda(double x);
	
	template <typename Pencil, typename Iterator, typename Segmentation>
	std::vector<Pencil> computeTangentialCover(Iterator itB, Iterator itE,
											   const Segmentation& s);


	template <typename Pencil, typename Iterator>
	std::vector<Pencil> theoreticalTangentsOnBoudin(Iterator itB, Iterator itE,double logFactor);


	template <typename Pencil, typename Iterator>
	std::vector<Pencil> orthogonalPlanesWithTangents(Iterator itb, Iterator ite);


	template <typename Pencil, typename Iterator>
	std::vector<Pencil> orthogonalPlanesWithNaiveTangents(Iterator itb, Iterator ite);

}

template <typename Pencil, typename Iterator>
std::vector<Pencil> TangentUtils::orthogonalPlanesWithTangents(Iterator itb, Iterator ite)
{
	typedef StandardDSS6Computer<Iterator,int,8> SegmentComputer;
	typedef SaturatedSegmentation<SegmentComputer> Segmentation;
    SegmentComputer algo;
    Segmentation s(itb, ite, algo);
    s.setMode("MostCentered++");
	typename Segmentation::SegmentComputerIterator i = s.begin();
	typename Segmentation::SegmentComputerIterator end = s.end();

	std::vector<Pencil> tangents = TangentUtils::computeTangentialCover<Pencil>(itb, ite, s);
	return tangents;
}

template <typename Pencil, typename Iterator>
std::vector<Pencil> TangentUtils::orthogonalPlanesWithNaiveTangents(Iterator itb, Iterator ite) {
	std::vector<Pencil> tangents;
	for (; itb != ite; ++itb) {
		auto nextIt = itb;
		++nextIt;
		Z3i::Point point = *itb;
		Z3i::Point tangent =  *nextIt - *itb;
		tangents.push_back({point, tangent});
	}
	return tangents;
}


/* Definition of the computation of tangential cover */
template <typename Pencil, typename Iterator, typename Segmentation>
std::vector<Pencil> TangentUtils::computeTangentialCover(Iterator itB, Iterator itE,
										   const Segmentation& s) {
	typedef typename Segmentation::SegmentComputerIterator SegmentComputerIterator;
	typedef typename Segmentation::SegmentComputer DSS;
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
					Vector3D v = *(--currentSegment.end()) - *currentSegment.begin();
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
		pencil.computeTangent(tangents);
		pencils.push_back(pencil);
	}
	return pencils;
}

template <typename Pencil, typename Iterator>
std::vector<Pencil> TangentUtils::theoreticalTangentsOnBoudin(Iterator itb, Iterator ite, double logFactor) {
	std::vector<Pencil> tangents;
	for (;itb != ite; ++itb) {
		typename Pencil::P point = *itb;
		double t = exp(point[2]/logFactor);
		double ztang = logFactor/t;
	    typename Pencil::Vector3d p1 = {1.+(double)point[0],1.+(double)point[1], ztang+(double)point[2]};
	    typename Pencil::Vector3d p2 = {-1.+(double)point[0], -1.+(double)point[1], -ztang+(double)point[2]};
		typename Pencil::Vector3d normal = (p1 - p2).getNormalized();
		tangents.push_back({point, normal});
	}
	return tangents;
}


#endif
