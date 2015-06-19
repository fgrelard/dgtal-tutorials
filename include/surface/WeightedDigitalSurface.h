#ifndef __WEIGHTED_DIGITAL_SURFACE__
#define __WEIGHTED_DIGITAL_SURFACE__

#include "DGtal/topology/DigitalSurface.h"
#include <vector>

template <typename Container, typename SurfacePoint>
class WeightedDigitalSurface : public DGtal::DigitalSurface<Container> {
public:
	typedef DGtal::DigitalSurface<Container> Base;
	typedef typename Base::Vertex Vertex;

	// -------- Standard services ----------- //
public:
	using Base::Base;
	WeightedDigitalSurface(const Container & container, const std::vector<SurfacePoint>& surfaceVoxels) : Base{container}, mySurfaceVoxels{surfaceVoxels} {}
	WeightedDigitalSurface(const Base& base, const std::vector<SurfacePoint>& surfaceVoxels) : Base(base), mySurfaceVoxels{surfaceVoxels} {}
	WeightedDigitalSurface(const WeightedDigitalSurface & other) : Base{other}, mySurfaceVoxels{other.mySurfaceVoxels} {}
	WeightedDigitalSurface(Container* containerPtr, const std::vector<SurfacePoint>& surfaceVoxels) : Base{containerPtr}, mySurfaceVoxels{surfaceVoxels}  {}

	// ---------- Undirected simple graph ---------- //
public:	
	template <typename OutputIterator>
	void writeNeighbors(OutputIterator & it, const Vertex& v) const override;

	template <typename OutputIterator, typename VertexPredicate>
	void writeNeighbors(OutputIterator & it, const Vertex& v, const VertexPredicate & pred) const override;

private:
	std::vector<SurfacePoint> mySurfaceVoxels;
};

template <typename Container, typename SurfacePoint>
template <typename OutputIterator>
void WeightedDigitalSurface<Container, SurfacePoint>::
writeNeighbors( OutputIterator & it, const Vertex & v) const {
	std::vector<Vertex> neighbors;
	OutputIterator tmp(neighbors);
	Base::writeNeighbors(tmp, v);
	for (auto itN = neighbors.begin(), itNE = neighbors.end(); itN!=itNE; ++itN) {
		//	*it++ = *itN;
		for (auto itSurfPoint = mySurfaceVoxels.begin(),
				 itSurfPointE = mySurfaceVoxels.end();
				 itSurfPoint != itSurfPointE; ++itSurfPoint) {
			SurfacePoint current = *itSurfPoint;
			if (current.contains(*itN)) {
				for (const auto& s : current.surfelsAtPoint()) {
						*it++ = s; 
				}
			}
		}
	}
}

template <typename Container, typename SurfacePoint>
template <typename OutputIterator, typename VertexPredicate>
void WeightedDigitalSurface<Container, SurfacePoint>::
writeNeighbors( OutputIterator & it, const Vertex & v, const VertexPredicate & pred) const {
	std::vector<Vertex> neighbors;
	OutputIterator tmp(neighbors);
	Base::writeNeighbors(tmp, v, pred);
	for (auto itN = neighbors.begin(), itNE = neighbors.end(); itN!=itNE; ++it) {
		*it++ = *itN;
		for (auto itSurfPoint = mySurfaceVoxels.begin(),
				 itSurfPointE = mySurfaceVoxels.end();
			 itSurfPoint != itSurfPointE; ++itSurfPoint) {
			SurfacePoint current = *itSurfPoint;
			if (current.contains(*itN)) {
				for (const auto& s : current.surfelsAtPoint()) { 
					if (pred(s)) *it++ = s; 
				}
			}
		}
	}
}

#endif
