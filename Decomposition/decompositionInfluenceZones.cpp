/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file dvcm-2d.cpp
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/04/15
 *
 * Computes the Voronoi Covariance Measure of a list of 2D digital
 * points. Displays the resulting normal vector and feature detection.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include <iterator>
#include <QtGui/qapplication.h>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "DGtal/math/linalg/EigenDecomposition.h"
#include "DGtal/geometry/volumes/distance/ExactPredicateLpSeparableMetric.h"
#include "geometry/VoronoiCovarianceMeasure.h"
#include "geometry/VoronoiCovarianceMeasureOnDigitalSurface.h"
#include "DGtal/io/colormaps/GradientColorMap.h"
#include "DGtal/io/boards/Board2D.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/shapes/ShapeFactory.h"
#include "DGtal/shapes/GaussDigitizer.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/topology/ImplicitDigitalSurface.h"
#include "DGtal/math/Histogram.h"
#include "DGtal/math/Statistic.h"
#include "DGtal/geometry/curves/StandardDSS6Computer.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"
#include "DGtal/graph/GraphVisitorRange.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "geometry/TangentUtils.h"
#include "geometry/SliceUtils.h"
#include "geometry/Pencil.h"
#include "geometry/MSTTangent.h"
#include "geometry/PointUtil.h"
#include "geometry/WeightedPoint.h"
#include "geometry/MedialAxis.h"
#include "geometry/ImageUtil.h"
#include "surface/SurfaceUtils.h"
#include "surface/Morphomaths.h"
#include "clustering/diana.hpp"
#include "geometry/WeightedPointCount.h"
#include "geometry/VCMUtil.h"
#include "geometry/WeightedPointCount.h"
#include "surface/SurfaceTraversal.h"
#include "Statistics.h"
#include "shapes/Ball.h"
///////////////////////////////////////////////////////////////////////////////
#include "ReverseHatPointFunction.h"

#include "graph/GraphEdge.h"
#include "graph/Concatenation.h"
#include "graph/LevelConcatenation.h"
#include "geometry/CurveAnalyzer.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;


Z3i::Point extractNearestNeighborInSetFromPoint(const Z3i::DigitalSet& aSet, const Z3i::RealPoint& aPoint) {
    double distanceMin = numeric_limits<double>::max();
    Z3i::Point toReturn;
    for (auto it = aSet.begin(), ite = aSet.end(); it != ite; ++it) {
        double distanceToPoint = sqrt(pow(((*it)[0] - aPoint[0]), 2) + pow(((*it)[1] - aPoint[1]), 2) + pow(((*it)[2] - aPoint[2]), 2));
        if (distanceToPoint < distanceMin || (distanceMin == distanceToPoint && aPoint > *it)) {
            distanceMin = distanceToPoint;
            toReturn = *it;
        }
    }
    return toReturn;
}


template <typename DTL2>
Z3i::Point findMaxDTInSet(const Z3i::DigitalSet& set, const DTL2 dt, const Z3i::Point& junctionPoint) {
    double maxDT = 0.0;
    Z3i::Point maxDTPoint;
    for (auto it = set.begin(), ite = set.end(); it != ite; ++it) {
        Z3i::Point point = *it;
        if (dt(point) > maxDT) {
            maxDT = dt(point);
            maxDTPoint = point;
        }
        else if (dt(point) == maxDT) {
            if (Z3i::l2Metric(junctionPoint, point) < Z3i::l2Metric(junctionPoint, maxDTPoint))
                maxDTPoint = point;
        }
    }
    return maxDTPoint;
}



template <typename DTL2>
Z3i::Point projectionPoint(const DTL2& dt,
                           const Z3i::Point& center,
                           const Z3i::RealVector& direction) {
    Z3i::Point proj = center;
    int scalar = 1;
    while (dt.domain().isInside(proj) && dt(proj) > 1) {
        scalar++;
        proj = center + direction*scalar;
    }
    return proj;
}

template <typename DTL2>
Z3i::RealVector signVector(const DTL2& dt,
                           const Z3i::Point& center,
                           const Z3i::RealVector& direction) {
    Z3i::Point projCenter = projectionPoint(dt, center, direction);
    Z3i::Point otherProjCenter = projectionPoint(dt,center, -direction);

    if (Z3i::l2Metric(projCenter, center) > Z3i::l2Metric(otherProjCenter, center))
        return -direction;
    return direction;
}




template <typename DTL2>
double deltaEdge(const DTL2& dt,
             const Z3i::DigitalSet& edge,
             const Z3i::DigitalSet& branch,
             const Z3i::DigitalSet& endPoints) {

    Z3i::Point b, end;
    for (const Z3i::Point& e : edge) {
        if (branch.find(e) != branch.end())
            b = e;
        if (endPoints.find(e) != endPoints.end())
            end = e;
    }
    double delta = 0;
    if (dt.domain().isInside(b) &&
        dt.domain().isInside(end) &&
        b != Z3i::Point() &&
        end != Z3i::Point())
        delta = dt(b) - dt(end);
    return delta;
}

template <typename DTL2>
int numberOfLocalMaximaDT(const DTL2& dt,
                          const Z3i::DigitalSet& aSet) {
    typedef MetricAdjacency<Z3i::Space, 3> MAdj;

    int localMax = 0;

    for (const Z3i::Point& p : aSet) {
        vector<Z3i::Point> neighbors;
        back_insert_iterator<vector<Z3i::Point>> inserter(neighbors);
        MAdj::writeNeighbors(inserter, p);
        double currentDTValue = dt(p);
        bool isLocalMaximum = true;
        for (const Z3i::Point& n: neighbors) {
            if (dt.domain().isInside(n) && dt(n) > currentDTValue)
                isLocalMaximum = false;
        }
        if (isLocalMaximum) localMax++;
    }
    return localMax;
}

unsigned int lengthEdge(const Z3i::DigitalSet& edge) {
    return edge.size();
}

template <typename DTL2>
double amountInformationLostRatio(const DTL2& dt,
                       const Z3i::DigitalSet& edge,
                       const Z3i::DigitalSet& branch,
                       const Z3i::DigitalSet& endPoints) {
    Z3i::Point b, end;
    for (const Z3i::Point& e : edge) {
        if (branch.find(e) != branch.end())
            b = e;
        if (endPoints.find(e) != endPoints.end())
            end = e;
    }
    double length = lengthEdge(edge);
    double delta = deltaEdge(dt, edge, branch, endPoints);
    double qratio = 0;
    if (dt.domain().isInside(b)) {
        double Rh = dt(b);
        double d = Z3i::l2Metric(b, end);
        qratio = length * (d - delta) / Rh;
    }
    return qratio;
}

void divideRegions() {

}

void ascribeRegions() {

}



///////////////////////////////////////////////////////////////////////////////
int main( int  argc, char**  argv )
{


    typedef Z3i::Space Space;
    typedef Z3i::Point Point;
    typedef Z3i::RealPoint RealPoint;
    typedef Z3i::RealVector RealVector;
    typedef HyperRectDomain<Space> Domain;
    typedef ExactPredicateLpSeparableMetric<Space, 2> Metric; // L2-metric
    typedef EigenDecomposition<3,double> LinearAlgebraTool;
    typedef LinearAlgebraTool::Matrix Matrix;
    typedef ImageSelector<Domain, unsigned char>::Type Image;
    typedef VoronoiCovarianceMeasure<Space,Metric> VCM;




    typedef WeightedPoint<Z3i::RealPoint> WeightedRealPoint;
    typedef WeightedPoint<Z3i::Point> WeightedPoint;
    typedef MetricAdjacency<Space, 3> MetricAdjacency;
    typedef WeightedPointCount<Point> WeightedPointCount;
    typedef functors::NotPointPredicate<Z3i::DigitalSet> NotPointPredicate;
    typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
    typedef DGtal::DistanceTransformation<Space, ThresholdedImage, Z3i::L2Metric> DTL2;
    typedef functors::BallConstantPointFunction<Point, double> KernelFunction;

    typedef Z3i::KSpace KSpace;
    typedef DGtal::functors::NotPointPredicate<ThresholdedImage> BackgroundPredicate;
    typedef ImplicitDigitalSurface< KSpace, BackgroundPredicate > DigitalSurfaceContainer;

    typedef VoronoiCovarianceMeasureOnDigitalSurface< DigitalSurfaceContainer, Metric, KernelFunction, DTL2> VCMOnSurface;
    typedef KSpace::Surfel Surfel;
    typedef VoronoiMap<Space, NotPointPredicate, Metric> VoronoiMap;
    typedef Eigen::MatrixXd MatrixXd;

    po::options_description general_opt("Allowed options are: ");
    general_opt.add_options()
        ("help,h", "display this message")
        ("input,i", po::value<std::string>(), "vol file (corresponding volume)")

        ("skeleton,s", po::value<std::string>(), "vol file (medial axis)")
        ("thresholdMin,m", po::value<int>()->default_value(0), "minimum threshold for binarization")
        ("thresholdMax,M", po::value<int>()->default_value(255), "maximum threshold for binarization")
        ;

    bool parseOK=true;
    po::variables_map vm;
    try{
        po::store(po::parse_command_line(argc, argv, general_opt), vm);
    } catch(const std::exception& ex){
        parseOK=false;
        trace.info()<< "Error checking program options: "<< ex.what()<< endl;
    }
    po::notify(vm);
    if( !parseOK || vm.count("help")||argc<=1)
    {
        std::cout << "Usage: " << argv[0] << " [input]\n"
                  << "Display volume file as a voxel set by using QGLviewer"<< endl
                  << general_opt << "\n";
        return 0;
    }
    if(!vm.count("input"))
    {
        trace.error() << " The file name was not defined" << endl;
        return 0;
    }


    string inputFilename = vm["input"].as<std::string>();
    string skeletonFilename = vm["skeleton"].as<std::string>();
    int thresholdMin = vm["thresholdMin"].as<int>();
    int thresholdMax = vm["thresholdMax"].as<int>();


    QApplication application(argc,argv);
    Viewer3D<> viewer;
    viewer.show();

    Image volume = VolReader<Image>::importVol(inputFilename);
    Z3i::Domain domainVolume = volume.domain();
    Z3i::DigitalSet setVolume(domainVolume);
    SetFromImage<Z3i::DigitalSet>::append<Image> (setVolume, volume,
                                                  thresholdMin-1, thresholdMax);
    Z3i::DigitalSet setSurface = SurfaceUtils::extractSurfaceVoxels(volume, thresholdMin, thresholdMax);


    Image skeleton = VolReader<Image>::importVol(skeletonFilename);
    Z3i::Domain domainSkeleton = skeleton.domain();
    Z3i::DigitalSet setSkeleton(domainSkeleton);

    SetFromImage<Z3i::DigitalSet>::append<Image>(setSkeleton, skeleton, thresholdMin-1, thresholdMax);

    Z3i::DigitalSet existingSkeleton = CurveAnalyzer::ensureConnexity(setSkeleton);
    typedef StandardDSS6Computer<vector<Point>::iterator,int,8> SegmentComputer;
    typedef GreedySegmentation<SegmentComputer> Segmentation;
        vector<Z3i::Point> endPointsV = CurveAnalyzer::findEndPoints(existingSkeleton);
    Z3i::Point p = (*endPointsV.begin());
    Z3i::DigitalSet branchingPoints(domainVolume);
    vector<Point> existingSkeletonOrdered = CurveAnalyzer::curveTraversalForGraphDecomposition(branchingPoints,
                                                                                               existingSkeleton,
                                                                                               p);


    //branchingPoints = CurveAnalyzer::detectCriticalPoints(existingSkeleton);
    vector<Z3i::Point> endPoints;
    for (const Z3i::Point& p : endPointsV) {
        if (branchingPoints.find(p) == branchingPoints.end())
            endPoints.push_back(p);
    }
    vector<Z3i::DigitalSet> edgeGraph = CurveAnalyzer::constructGraph(existingSkeletonOrdered, branchingPoints);
    Z3i::DigitalSet endPointSet(domainSkeleton);
    endPointSet.insert(endPoints.begin(), endPoints.end());
    vector<GraphEdge*> hierarchicalGraph = CurveAnalyzer::hierarchicalDecomposition(edgeGraph, endPoints, branchingPoints);

    Metric l2;
    typedef DGtal::functors::IntervalForegroundPredicate<Image> Binarizer;
    Binarizer binarizer(volume, thresholdMin-1, thresholdMax);
    DTL2 dt(&volume.domain(), &binarizer, &Z3i::l2Metric);

    vector<Z3i::DigitalSet> zonesOfInfluence;
    for (const Z3i::Point& b : branchingPoints) {
            Z3i::DigitalSet zoneOfInfluence(setVolume.domain());
            double radius = dt(b);
            Ball<Z3i::Point> ball(b, radius);
            vector<Z3i::Point> pointsInBall = ball.intersection(setVolume);
            zoneOfInfluence.insert(pointsInBall.begin(), pointsInBall.end());
            for (const Z3i::Point& p : pointsInBall) {
                    double radius = dt(p);
                    Ball<Z3i::Point> ball(p, radius);
                    vector<Z3i::Point> pointsInBall = ball.intersection(setVolume);
                    zoneOfInfluence.insert(pointsInBall.begin(), pointsInBall.end());
            }
            zonesOfInfluence.push_back(zoneOfInfluence);
    }



    for (const Z3i::DigitalSet& zone : zonesOfInfluence) {
        viewer << CustomColors3D(Color::Red, Color::Red) << zone;
    }

    viewer << CustomColors3D(Color::Red, Color::Red) << existingSkeleton;
    viewer << CustomColors3D(Color(210,210,210,20), Color(210,210,210,20)) << setVolume;

    //second pass



    viewer << Viewer3D<>::updateDisplay;
    application.exec();
    return 0;
}
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
