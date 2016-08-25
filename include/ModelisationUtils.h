#ifndef MODELISATION_UTILS_H
#define MODELISATION_UTILS_H

#include <vector>
#include <set>
#include <Eigen/Dense>
#include <Eigen/Geometry>

template <typename Point>
void add(const Point& p, std::vector<Point> & v);


/**
 * This function finds if the 6 neighbour of a given point p lie on the same orientation
 * If there are more than 2 direct neighbours, then we don't have the same orientation
 * that is to say we dont have a curve
 * This function allows to determine points which can be removed
 **/

template <typename Point>
bool sameDirection6Neighbour(const std::vector<Point> & v, Point & p);

/**
 * Allows to obtain a 26 connected curve
 **/
template <typename Point>
void construct26ConnectedCurve(std::vector<Point> & vectorPoints);

template <typename Point>
void drawCircle(std::vector<Point> & vectorPoints, float radius, float cx, float cy, float z, float increment);

void drawDisk(Eigen::Matrix<double, Eigen::Dynamic, 4> & m, double radius, double cx, double cy, double cz, long int &row, float increment);

template <typename Point>
void drawDisk(std::vector<Point>& v, float radius, float cx, float cy, float cz, float increment);

template <typename Point>
void drawCone(std::vector<Point>& v, int length, float radius, float increment);

template <typename Point>
void drawCylinder(std::vector<Point>& v, int length, float radius, float increment);

void drawCylinder(Eigen::Matrix<double, Eigen::Dynamic, 4> & m, int length, int radius, float cz, float rotationAngle, float increment);


/*
 * Radius  Winding : radius of the loop
 * RadiusSpiral: radius of the volume
 */
template <typename Point>
void createHelixCurve(std::vector<Point> & vectorPoints, int range, int radiusWinding, int radiusSpiral, int pitch, float increment);


/*
 * This function does not allow to set the spiral radius
 */
template <typename Point>
void createHelixCurve( std::vector<Point> &point, int range, int radius, int pitch, float increment);

template <typename Point>
void createStraightLine(std::vector<Point> & points, int range);

template <typename Point>
void createSyntheticAirwayTree(std::vector<Point> & v, int branchNumber, int lengthTrachea, int z, float rotationAngle, Point firstPoint, float increment);

template <typename Point>
void createLogarithmicCurve(std::vector<Point> & curve, int range, float increment);

template <typename CurvePoint, typename Point>
void createVolumeFromCurve(const std::vector<Point> & curve, std::set<Point> & volume, int ballRadius);

template <typename CurvePoint, typename Point>
void createRotatedVolumeFromCurve(const std::vector<CurvePoint> & curve, std::set<Point> & volume, int ballRadius, double angle, const Eigen::Vector3d& vector = Eigen::Vector3d::UnitX());

template <typename DigitalSet>
DigitalSet addNoise(const DigitalSet& set);

template <typename DigitalSet>
void create2DCurve(DigitalSet & digitalSet);

template <typename DigitalSet, typename Board>
void create2DNaiveTangentsForVisu(const DigitalSet & points, Board& board);

#include "ModelisationUtils.ih"

#endif
