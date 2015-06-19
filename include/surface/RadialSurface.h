#ifndef __RADIAL_DIGITAL_SURFACE__
#define __RADIAL_DIGITAL_SURFACE__

#include "DGtal/topology/DigitalSurface.h"
#include <vector>

template <typename Container>
class RadialSurface : public DGtal::DigitalSurface<Container> {
public:
	typedef DGtal::DigitalSurface<Container> Base;
	typedef typename Base::Vertex Vertex;
	typedef typename Base::KSpace KSpace;

	// -------- Standard services ----------- //
public:
	using Base::Base;
	RadialSurface(const Container & container, int aForbiddenDirection) : Base{container}, myForbiddenDirection{aForbiddenDirection} {}
	RadialSurface(const Base& base, int aForbiddenDirection) : Base{base}, myForbiddenDirection{aForbiddenDirection} {}
	RadialSurface(const RadialSurface & other) : Base{other}, myForbiddenDirection{other.myForbiddenDirection} {}

	// ---------- Undirected simple graph ---------- //
public:	
	template <typename OutputIterator>
	void writeNeighbors(OutputIterator & it, const Vertex& v) const override;

private:
    int myForbiddenDirection;
};

template <typename Container>
template <typename OutputIterator>
void RadialSurface<Container>::
writeNeighbors( OutputIterator & it, const Vertex & v) const {
	Vertex s;
	auto aTracker = Base::container().newTracker(v);
	aTracker->move(v);
	for ( typename KSpace::DirIterator q = Base::container().space().sDirs( v );
		  q != 0; ++q )
    {
		if ( aTracker->adjacent( s, *q, true ) )
			*it++ = s;
		if ( aTracker->adjacent( s, *q, false ) ) {
			if (*q == myForbiddenDirection) continue;
			*it++ = s;
		}
    }
}

#endif
