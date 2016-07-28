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
                WatershedPoint() {}
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

                void setValue(double aValue) { myValue = aValue; }
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
                bool operator()(T* lhs, T* rhs) {
                        return *lhs < *rhs;
                }
        };



public:
        template <typename Container>
        Watershed(const Container& container, double aEpsilon);

        void compute();
        void pixelsAtSameAltitude(std::queue<Point>& fifo, Point& current, double currentAltitude);
        void extendBasins(std::queue<Point>& fifo);
        void detectMinima(std::queue<Point>& fifo, int& label, double altitude);


        std::map<Point, WatershedInformation*> getWatershed();
        int getBins();
private:
        std::map<Point, WatershedInformation*> myImageWatershed;
        std::map<WatershedInformation*, Point, ptr_less<WatershedInformation>> sortedMap;
        double myEpsilon;
        Point myFictitious;
};


template <typename Point>
template <typename Container>
Watershed<Point>::Watershed(const Container& container, double aEpsilon) {
        myEpsilon = aEpsilon;
        for (const auto& p : container) {
                WatershedInformation* wi = new WatershedInformation(p.second);
                myImageWatershed[p.first] = wi;
                sortedMap[wi] = p.first;
                if (p.first < myFictitious)
                        myFictitious = 2*p.first;
        }
}





template <typename Point>
void Watershed<Point>::compute() {
  int label = UNLABELLED;
  std::queue<Point> fifo;

  Point current = myImageWatershed.begin()->first;
  DGtal::trace.beginBlock("Watershed");
  double altitude = 0, previousAltitude = -1;
  while ( altitude != previousAltitude ) {
          previousAltitude = altitude;
          altitude = myImageWatershed[current]->getValue();
          trace.info() << altitude << endl;
          pixelsAtSameAltitude(fifo, current, altitude);
          extendBasins(fifo);
          detectMinima(fifo, label, altitude);
  }
  DGtal::trace.endBlock();
//  closeWatershed(t_pixels);
}



template <typename Point>
void Watershed<Point>::pixelsAtSameAltitude(std::queue<Point>& fifo, Point& current, double currentAltitude) {
        DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;
        for (const std::pair<WatershedInformation*, Point>& pairWatershed : sortedMap) {
                Point p = pairWatershed.second;
                WatershedInformation* wi = myImageWatershed[p];
                if (wi->getValue() > currentAltitude + myEpsilon) break;
                wi->setLabel(MASK);
                std::vector<Point> neighbors;
                std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
                adj.writeNeighbors(inserter, p);
                for (const Point& n : neighbors) {
                        if (myImageWatershed.find(n) != myImageWatershed.end()) {
                                WatershedInformation* wn = myImageWatershed.at(n);
                                if (wn->getLabel() > 0 || wn->getLabel() == WATERSHED) {
                                        wi->setDistance(1);
                                        fifo.push(p);
                                        break;
                                }
                        }
                }
                current = p;
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

                WatershedInformation* wp = myImageWatershed[p];
                for (const Point& n : neighbors) {
                        if (myImageWatershed.find(n) != myImageWatershed.end()) {
                                WatershedInformation* wn = myImageWatershed.at(n);
                                if (wn->getDistance() <= d_cur &&
                                    (wn->getLabel() == WATERSHED || wn->getLabel() > 0)) {
                                        if (wn->getLabel() > 0) {
                                                if (wp->getLabel() == MASK)
                                                        wp->setLabel(wn->getLabel());
                                                else if (wp->getLabel() != wn->getLabel())
                                                        wp->setLabel(WATERSHED);
                                        }
                                        else if (wp->getLabel() == MASK) {
                                                wp->setLabel(WATERSHED);
                                        }
                                }
                                else if (wn->getLabel() == MASK && wn->getDistance() == 0) {
                                        wn->setDistance(d_cur+1);
                                        fifo.push(n);
                                }
                        }
                }

        } // End : while ( true )
}


template <typename Point>
void Watershed<Point>::detectMinima(std::queue<Point>& fifo, int& label, double altitude) {
        DGtal::MetricAdjacency<DGtal::Z3i::Space, 3> adj;
        for (std::pair<Point, WatershedInformation*> pairWatershed : myImageWatershed) {
                Point p = pairWatershed.first;
                WatershedInformation* wp = pairWatershed.second;
                if (wp->getValue() >= altitude && wp->getValue() <= altitude + myEpsilon) {
                        wp->setDistance(0);
                        if (wp->getLabel() == MASK) {
                                label++;
                                fifo.push(p);
                                wp->setLabel(label);
                                trace.info() << p << " " << wp->getLabel() << endl;
                                while (!fifo.empty()) {
                                        Point q = fifo.front();
                                        fifo.pop();
                                        // Labelling p by inspecting neighbours
                                        std::vector<Point> neighbors;
                                        std::back_insert_iterator<std::vector<Point> > inserter(neighbors);
                                        adj.writeNeighbors(inserter, q);
                                        for (const Point& n : neighbors) {
                                                if (myImageWatershed.find(n) != myImageWatershed.end()) {
                                                        WatershedInformation* wn = myImageWatershed.at(n);
                                                        if (wn->getLabel() == MASK) {
                                                                fifo.push(n);
                                                                wn->setLabel(label);
                                                        }
                                                }
                                        }

                                }
                        }
                }
        }
}

template <typename Point>
std::map<Point, typename Watershed<Point>::WatershedInformation*> Watershed<Point>::getWatershed() {
        return myImageWatershed;
}

template <typename Point>
int Watershed<Point>::getBins() {
        return (*max_element(myImageWatershed.begin(), myImageWatershed.end(), [](const std::pair<Point, WatershedInformation*>& one,
                                                                                  const std::pair<Point, WatershedInformation*>& two) {
                                     return one.second->getLabel() < two.second->getLabel();
                             })).second->getLabel();
}



#endif
