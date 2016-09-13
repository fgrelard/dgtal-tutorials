#ifndef DIGITALPLANE_H
#define DIGITALPLANE_H

#include "DGtal/geometry/surfaces/ParallelStrip.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/kernel/sets/DigitalSetSelector.h"
#include "DGtal/topology/Object.h"
#include "DGtal/topology/DigitalTopology.h"
#include "DGtal/topology/MetricAdjacency.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
template <typename TSpace>
class DigitalPlane
{
public :
        typedef DGtal::ParallelStrip<TSpace> PlaneEquation;
        typedef typename TSpace::Point Point;
        typedef typename TSpace::RealVector Vector;
        typedef typename DGtal::HyperRectDomain< TSpace > Domain;
        typedef typename DGtal::DigitalSetSelector< Domain, DGtal::BIG_DS+DGtal::HIGH_BEL_DS >::Type DigitalSet;
        typedef DGtal::MetricAdjacency< TSpace, 1 > Adj4;
        typedef DGtal::MetricAdjacency< TSpace, 2 > Adj8;
        typedef DGtal::MetricAdjacency< TSpace, 3> Adj26;
        typedef DGtal::DigitalTopology<Adj26, Adj4> Topology;
        typedef DGtal::ExactPredicateLpSeparableMetric<TSpace, 2> L2Metric;

public:
        DigitalPlane(const Point& aPoint, const Vector& aNormal, int aConnexity = 26) : myPoint(aPoint) {
                double omega, d;
                for (int i = 0; i < TSpace::dimension; i++)
                        d += aNormal[i] * aPoint[i];
                if (aConnexity == 26 || aConnexity == 8) {
                        omega = *std::max_element(aNormal.begin(), aNormal.end());
                }
                else // if (aConnexity == 6 || aConnexity == 4)
                {
                        for (auto it = aNormal.begin(), ite = aNormal.end(); it != ite; ++it)
                                omega += *it;
                }
                myPlaneEquation = PlaneEquation(d, aNormal, omega);
        }


        DigitalSet intersectionWithSet(const DigitalSet& pointsV) {
                typedef typename DigitalSet::Point Value;
                DigitalSet points(pointsV);
                points.clear();

                for (const Value& value : pointsV) {
                        if (contains(value))
                                points.insert(value);
                }
                return points;
        }

        DigitalSet intersectionWithSetOneCC(const DigitalSet& pointsV) {
                typedef DGtal::Object<Topology, DigitalSet> ObjectType;
                L2Metric l2Metric;
                Adj26 adj26;
                Adj4 adj4;
                Topology dt26_6(adj26, adj4, DGtal::JORDAN_DT);
                DigitalSet intersection = intersectionWithSet(pointsV);
                ObjectType objectIntersection(dt26_6, intersection);
                std::vector<ObjectType> objects;
                std::back_insert_iterator<std::vector<ObjectType> > inserter(objects);
                unsigned int nbConnectedComponents = objectIntersection.writeComponents(inserter);
                DigitalSet connectedComponent = intersection;
                double min = std::numeric_limits<double>::max();
                if (nbConnectedComponents > 1) {
                        for (auto it = objects.begin(), ite = objects.end(); it != ite; ++it) {
                                double sum = 0;
                                DigitalSet ccSet = it->pointSet();
                                for (auto it = ccSet.begin(), ite = ccSet.end(); it != ite; ++it) {
                                        sum += l2Metric(*it, myPoint);
                                }
                                sum /= ccSet.size();
                                if (sum < min) {
                                        min = sum;
                                        connectedComponent = ccSet;
                                }
                        }
                }
                return connectedComponent;
        }

        bool contains(const Point& aPoint) {
                double valueToCheck;
                double omega = myPlaneEquation.nu();
                double d = myPlaneEquation.mu();
                Vector normal = myPlaneEquation.normal();
                for (int i = 0; i < TSpace::dimension; i++)
                        valueToCheck += aPoint[i] * normal[i];
                if (valueToCheck >= d && valueToCheck < d + omega) {
                        return true;
                }
                return false;
        }

        bool isPointAbove(const Point& aPoint) {
                double d = myPlaneEquation.mu();
                Vector normal = myPlaneEquation.normal();

                double valueToCheck;
                for (int i = 0; i < TSpace::dimension; i++)
                        valueToCheck += aPoint[i] * normal[i];
                if (valueToCheck >= d)
                        return true;
                return false;
        }

private:
        PlaneEquation myPlaneEquation;
        Point myPoint;
};

#endif
