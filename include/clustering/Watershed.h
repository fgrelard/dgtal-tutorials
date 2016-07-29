#ifndef WATERSHED_H
#define WATERSHED_H

#include "geometry/WeightedPointCount.h"
#include <vector>
#include <queue>
#include <map>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/topology/MetricAdjacency.h"
#include <boost/bimap.hpp>

#define MASK -3
#define WATERSHED -2
#define UNLABELLED -1
#define FIRST_LABEL 0


template <typename Point>
class Watershed {


        class WatershedPoint : public WeightedPointCount<Point> {
                typedef WeightedPointCount<Point> Base;
        public:
                using Base::Base;
        public:
                WatershedPoint() : myDistance(-1) {}
                WatershedPoint(const Point& aPoint, double aWeight, int aCount, double aDistance) : Base(aPoint, aWeight, aCount), myDistance(aDistance) {}
                WatershedPoint(const WatershedPoint& other) : Base(other), myDistance(other.myDistance) {}

                friend bool operator<(const WatershedPoint& it, const WatershedPoint& other) {
                        return (it.myWeight < other.myWeight);
                }
        public:
                double myDistance;
        };


        class WatershedInformation {
        public:
                WatershedInformation() : myLabel(UNLABELLED), myDistance(0.0), myValue(0.0) {}
                WatershedInformation(double aValue) : myLabel(UNLABELLED), myDistance(0.0), myValue(aValue) {}

                void setLabel(int aLabel) { myLabel = aLabel; }
                int getLabel() const { return myLabel; }

                double getValue() const { return myValue; }

                void setDistance(double aDistance) { myDistance = aDistance; }
                double getDistance() const { return myDistance; }

                friend bool operator<(WatershedInformation it, WatershedInformation other) {
                        return (it.getValue() < other.getValue());
                }
        private:
                int myLabel;
                double myDistance;
                double myValue;
        };

        template<class T> struct ptr_less {
                bool operator()(const T* lhs, const T* rhs) {
                        return *lhs < *rhs;
                }
        };



public:
        template <typename Container>
        Watershed(const Container& container, double aEpsilon);

        void compute();
        void pixelsAtSameAltitude(std::queue<Point>& fifo, int& currentI, double currentAltitude);
        void extendBasins(std::queue<Point>& fifo);
        void detectMinima(std::queue<Point>& fifo, int& label, double altitude);

        WatershedPoint * at(const Point& p);

        std::vector<WatershedPoint*> getWatershed();
        int getBins();
private:
        // std::map<Point, WatershedInformation*> myImageWatershed;
        // std::map<WatershedInformation*, Point, ptr_less<WatershedInformation>> sortedMap;
        std::vector<WatershedPoint*> myImageWatershed;
        double myEpsilon;
        Point myFictitious;
};


template <typename Point>
template <typename Container>
Watershed<Point>::Watershed(const Container& container, double aEpsilon) {
        myEpsilon = aEpsilon;
        for (const auto& p : container) {
                WatershedPoint* wp = new WatershedPoint(p.first, p.second, UNLABELLED, 0.0);
                myImageWatershed.push_back(wp);
                // sortedMap[wi] = p.first;
                // if (p.first < myFictitious)
                //         myFictitious = 2*p.first;
        }
        sort(myImageWatershed.begin(), myImageWatershed.end(), ptr_less<WatershedPoint>());
}


template <typename Point>
typename Watershed<Point>::WatershedPoint* Watershed<Point>::at(const Point& p) {
        auto iterator = find_if(myImageWatershed.begin(), myImageWatershed.end(), [&](WatershedPoint* one) {
                        return one->myPoint == p;
                });
        if (iterator != myImageWatershed.end())
                return *iterator;
        else
                return new WatershedPoint();
}


template <typename Point>
void Watershed<Point>::compute() {
  int label = UNLABELLED;
  std::queue<Point> fifo;

  int currentI = 0;
  DGtal::trace.beginBlock("Watershed");
  double altitude = 0, previousAltitude = -1;
  while ( currentI < myImageWatershed.size() ) {
          previousAltitude = altitude;
          altitude = myImageWatershed[currentI]->myWeight;
          trace.info() << altitude << endl;
          pixelsAtSameAltitude(fifo, currentI, altitude);
          extendBasins(fifo);
          detectMinima(fifo, label, altitude);
  }
  DGtal::trace.endBlock();
//  closeWatershed(t_pixels);
}



template <typename Point>
void Watershed<Point>::pixelsAtSameAltitude(std::queue<Point>& fifo, int& currentI, double currentAltitude) {
        DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;
        while (currentI < myImageWatershed.size() &&
               myImageWatershed[currentI]->myWeight <= currentAltitude + myEpsilon) {
                WatershedPoint* wi = myImageWatershed[currentI];
                Point p = wi->myPoint;
                wi->myCount =MASK;
                std::vector<Point> neighbors;
                std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
                adj.writeNeighbors(inserter, p);
                for (const Point& n : neighbors) {
                        WatershedPoint* wn = at(n);
                        if (wn->myDistance != -1) {
                                if (wn->myCount > 0 || wn->myCount == WATERSHED) {
                                        wi->myDistance= 1;
                                        fifo.push(p);
                                        break;
                                }
                        }
                }
                currentI++;
        }
}

template <typename Point>
void Watershed<Point>::extendBasins(std::queue<Point>& fifo) {
  // The first pixels treated are the closest (d=1), then d=2...
        int d_cur = 1;
        fifo.push(myFictitious);
        DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;

        while ( !fifo.empty() ) {
                Point p = fifo.front();
                fifo.pop();
                if (p == myFictitious)  {
                        if ( fifo.empty() ) // this altitude is processed
                                break;
                        else {
                                fifo.push(myFictitious);
                                d_cur++;
                                p=fifo.front();
                                fifo.pop();
                        }
                }
                // Labelling p by inspecting neighbours
                std::vector<Point> neighbors;
                std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
                adj.writeNeighbors(inserter, p);

                WatershedPoint* wp = at(p);
                for (const Point& n : neighbors) {
                        WatershedPoint* wn = at(n);
                        if (wn->myDistance != -1) {
                                if (wn->myDistance <= d_cur &&
                                    (wn->myCount == WATERSHED || wn->myCount > 0)) {
                                        if (wn->myCount > 0) {
                                                if (wp->myCount == MASK)
                                                        wp->myCount =wn->myCount;
                                                else if (wp->myCount != wn->myCount)
                                                        wp->myCount =WATERSHED;
                                        }
                                        else if (wp->myCount == MASK) {
                                                wp->myCount =WATERSHED;
                                        }
                                }
                                else if (wn->myCount == MASK && wn->myDistance == 0) {
                                        wn->myDistance = d_cur+1;
                                        fifo.push(n);
                                }
                        }
                }

        } // End : while ( true )
}


template <typename Point>
void Watershed<Point>::detectMinima(std::queue<Point>& fifo, int& label, double altitude) {
        DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;
        for (WatershedPoint * wp : myImageWatershed) {
                Point p = wp->myPoint;
                if (wp->myWeight >= altitude && wp->myWeight <= altitude + myEpsilon) {
                        wp->myDistance = 0;
                        if (wp->myCount == MASK) {
                                label++;
                                fifo.push(p);
                                wp->myCount =label;
                                while (!fifo.empty()) {
                                        Point q = fifo.front();
                                        fifo.pop();
                                        // Labelling p by inspecting neighbours
                                        std::vector<Point> neighbors;
                                        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
                                        adj.writeNeighbors(inserter, q);
                                        for (const Point& n : neighbors) {
                                                WatershedPoint * wn = at(n);
                                                if (wn->myDistance != -1) {
                                                        if (wn->myCount == MASK) {
                                                                fifo.push(n);
                                                                wn->myCount = label;
                                                        }
                                                }
                                        }

                                }
                        }
                }
        }
}

template <typename Point>
std::vector<typename Watershed<Point>::WatershedPoint*> Watershed<Point>::getWatershed() {
        return myImageWatershed;
}

template <typename Point>
int Watershed<Point>::getBins() {
        return (*max_element(myImageWatershed.begin(), myImageWatershed.end(), [](WatershedPoint* one,
                                                                                  WatershedPoint* two) {
                                     return one->myCount < two->myCount;
                             }))->myCount;
}



#endif
