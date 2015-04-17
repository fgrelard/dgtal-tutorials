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
//LICENSE-END
/**
 * @file testSegmentation.cpp
 * @ingroup Tests
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 *
 * @date 2011/07/22
 *
 * This file is part of the DGtal library
 */

/**
 * Description of testSegmentation <p>
 * Aim: simple test of \ref GreedySegmentation and \ref SaturatedSegmentation
 */

#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <iterator>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/sets/DigitalSetInserter.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include <QtGui/qapplication.h>
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/base/Exceptions.h"

#include "DGtal/io/boards/Board2D.h"
#include "DGtal/io/Color.h"
#include "DGtal/geometry/curves/GridCurve.h"
#include "DGtal/base/Circulator.h"
#include "DGtal/topology/KhalimskySpaceND.h"
#include "DGtal/topology/SurfelAdjacency.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"
#include "DGtal/topology/DigitalSurface2DSlice.h"
#include "DGtal/topology/helpers/Surfaces.h"

#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/FreemanChain.h"
#include "DGtal/geometry/curves/SaturatedSegmentation.h"
#include "DGtal/geometry/curves/GreedySegmentation.h"


#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <Eigen/Dense>

#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/graph/BreadthFirstVisitor.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/LightImplicitDigitalSurface.h"
#include "DGtal/kernel/BasicPointPredicates.h"
#include "RosenProffittLengthEstimator.h"

// Local includes
#include "SurfacePoint.h"
#include "WeightedDigitalSurface.h"
#include "../CreateCurve/Distance.h"
#include "../Tangent3D/SliceUtils.h"
#include "../Tangent3D/Pencil.h"
#include "../Tangent3D/MSTTangent.h"
#include "Statistics.h"
#include "Ellipse.h"

using namespace DGtal;
using namespace std;
using namespace Z3i;
namespace po = boost::program_options;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class GreedySegmentation
///////////////////////////////////////////////////////////////////////////////


template <typename BreadthFirstVisitor, typename Surfel, typename Viewer, typename DigitalSurface>
void shortestPointPath(const DigitalSurface & digSurf, const Surfel& bel, const Surfel& end, Viewer& viewer) {
	typedef RosenProffittLengthEstimator<vector<Surfel>> LengthEstimator;
	map<Surfel, Surfel> aMapPointPrevious;
	BreadthFirstVisitor pvisitor(digSurf, bel);
	while (!pvisitor.finished()) {
		typename BreadthFirstVisitor::Node node = pvisitor.current();
		vector<Surfel> neighbors;
		back_insert_iterator<vector<Surfel>> iter(neighbors);
		pvisitor.graph().writeNeighbors(iter, node.first);
		for (auto it = neighbors.begin(), itE = neighbors.end(); it != itE; ++it) {
			
			auto itMapExisting = aMapPointPrevious.find(*it);
			if (itMapExisting == aMapPointPrevious.end())
			{
				aMapPointPrevious[*it] = node.first;
				if (*it == end) break;
			}
			else {
				//Compute Rosen Proffitt distance with previous map
				Surfel tmp = *it;
				vector<Surfel> path1;
				while (tmp != bel) {
					path1.push_back(tmp);
					tmp = aMapPointPrevious[tmp];
				}
				
				Surfel tmp2 = node.first;
				vector<Surfel> path2;
				path2.push_back(*it);
				while (tmp2 != bel) {
					path2.push_back(tmp2);
					tmp2 = aMapPointPrevious[tmp2];
				}
				LengthEstimator lengthPath;
				double d1 = lengthPath.eval(path1.begin(), path1.end());
				double d2 = lengthPath.eval(path2.begin(), path2.end());
				if (d1 > d2) {
					if (node.first != bel) {
					    aMapPointPrevious[*it] = node.first;
					}
				}
			}
		}
		pvisitor.expand();
				
	}
	std::cout << "visu" << std::endl;
	visualizePath(bel, end, aMapPointPrevious, viewer, Color::Green);

}

template <typename Surfel>
void visualizePath(const Surfel& start, const Surfel& end, map<Surfel, Surfel>& previous, Viewer3D<>& viewer, const Color& color) {
	Surfel tmp = end;
	while (tmp != start) {
		viewer << CustomColors3D(color, color) << tmp;
	    tmp = previous[tmp];
	}

}

template <typename Vector, typename Point>
Vector computePlaneNormal(const std::vector<Point> & points) {
	typedef Eigen::Matrix<double, Eigen::Dynamic, 3> MatrixXi;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXi;
	unsigned int size = points.size();
	MatrixXi A(size, 3);
	VectorXi b = VectorXi::Zero(size, 1);
	
	for (int i = 0; i < size; i++) {
		A(i, 0) = (double)points[i][0]*1.0;
		A(i, 1) = (double)points[i][1]*1.0;
		A(i, 2) = 1.0;
		b(i, 0) = (double)points[i][2]*1.0;
	}
	Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
	Vector normal;
	normal[0] = x(0, 0);
	normal[1] = x(1, 0);
	normal[2] = -1;
	return normal.getNormalized();
}


template <typename Point>
bool areAlmostSimilar(const Point& point, const Point& other) {
	typename Point::Scalar otherx = other[0];
	typename Point::Scalar othery = other[1];
	typename Point::Scalar otherz = other[2];

	typename Point::Scalar pointx = point[0];
	typename Point::Scalar pointy = point[1];
	typename Point::Scalar pointz = point[2];
	
	bool sameX = pointx == otherx || pointx == otherx + 1 || pointx == otherx-1;
	bool sameY = pointy == othery || pointy == othery + 1 || pointy == othery-1;
	bool sameZ = pointz == otherz || pointz == otherz + 1 || pointz == otherz-1;
	return sameX && sameY && sameZ;
}

template <typename Point>
bool areAlmostSimilar(const set<Point>& s1, const set<Point>& s2) {
	int cpt = 0;
	for (const Point& p1 : s1) {
		for (const Point& p2 : s2) {
			if (areAlmostSimilar(p1, p2)) {
				cpt++;
				break;
			}
		}
	}
	return cpt == s1.size();
}
template <typename Surfel, typename Point>
bool checkSymmetry(const map<Surfel, Point>& surfelToPoint, const set<Surfel>& path1, const set<Surfel>& path2, const Point & center) {
	int foundOneSymmetry = 0;
	for (auto it = path1.begin(), itE = path1.end(); it != itE; ++it) {
		Point currentPoint = surfelToPoint.at(*it);
		Point vectorToCenter = center - currentPoint;
		Point symmetryCurrent = center + vectorToCenter;
		for (auto itS = path2.begin(), itSE = path2.end(); itS != itSE; ++itS) {
			Point putativeSymmetric = surfelToPoint.at(*itS);
			if (areAlmostSimilar(putativeSymmetric, symmetryCurrent) && !areAlmostSimilar(center, symmetryCurrent) && !areAlmostSimilar(center, currentPoint)) {
				foundOneSymmetry++;
				break;
			}
		}
	}
	return (foundOneSymmetry == path1.size());
}

template <typename Point>
bool multipleDirection(const vector<Point>& aVector, int& cpt) {
	bool xOrientation = false, yOrientation = false, zOrientation = false;
    int posx = 0, posy = 0, posz = 0, negx = 0, negy = 0, negz = 0;
	for (auto it = aVector.begin(), itE = aVector.end(); it != itE; ++it) {
		Point diff = *it;
		if (diff[0] > 0) {
			posx++;
		}
		if (diff[0] < 0) {
			negx++;
		}
		if (diff[1] > 0) {
		    posy++;
		}
		if (diff[1] < 0) {
			negy++;
		}
		if (diff[2] > 0) {
			posz++;
		}
		if (diff[2] < 0) {
			negz++;
		}
	}
	xOrientation = (posx > 0 && negx > 0);
	yOrientation = (posy > 0 && negy > 0);
	zOrientation = (posz > 0 && negz > 0);
	if (xOrientation && yOrientation)
		cpt = posz > 0 ? posz : negz;
	else if (xOrientation && zOrientation)
		cpt = posy > 0 ? posy : negy;
	else if (zOrientation && yOrientation)
		cpt = posx > 0 ? posx : negx;
	return (xOrientation && yOrientation) || (xOrientation && zOrientation) || (zOrientation && yOrientation);
}

template <typename SurfacePoint, typename KSpace, typename Surfel, typename Point>
vector<SurfacePoint> computeSurfelWeight(const KSpace & ks, const DigitalSet & set3d, const set<Surfel>& surfelSet, map<Surfel, Point> & surfaceVoxelSet) {
	vector<SurfacePoint> weightedSurfaceVector;
	for (auto it = set3d.begin(), itE = set3d.end(); it != itE; ++it) {
		vector<Surfel> aSurfelV;
		SCell current = ks.sSpel(*it);
		int number = 0;
		for (int i = 0; i < 3; i++) {
			auto itSurfel = surfelSet.find(ks.sIncident(current, i, true));
			if (itSurfel != surfelSet.end()) {
				number++;
				aSurfelV.push_back(*itSurfel);
				surfaceVoxelSet[*itSurfel] = *it; 
			}
		}

		for (int i = 0; i < 3; i++) {
			auto itSurfel = surfelSet.find(ks.sIncident(current, i, false));
			if (itSurfel != surfelSet.end()) {
				number++;
				aSurfelV.push_back(*itSurfel);
				surfaceVoxelSet[*itSurfel] = *it;
			}
		}
		weightedSurfaceVector.push_back({*it, number, aSurfelV});
		// if (number > 0)
		// 	surfaceVoxelSet.insert(ks.sCoords(current));
	}
	return weightedSurfaceVector;
}

template <typename Point>
vector<Point> nearestPointsFromCenter(const vector<Point> & points, const Point & center, double radius) {
	vector<Point> nearestPointV;
	for (const Point& point : points) {
		if (euclideanDistance(point, center) <= radius) {
			nearestPointV.push_back(point);
		}
	}
	return nearestPointV;
}

template <typename Surfel, typename Point>
bool areNeighbours(const set<Surfel>& path1, const set<Surfel>& path2, const map<Surfel, Point>& surfelToPoint) {
	int cpt = 0;
	for (const auto & surfelInPOne : path1) {
		Point pp1 = surfelToPoint.at(surfelInPOne);
		for (const auto & surfelInPTwo : path2) {
			Point pp2 = surfelToPoint.at(surfelInPTwo);
			if (areAlmostSimilar(pp1, pp2))
				cpt++;
		}
	}
	trace.info() << cpt <<endl;
	return cpt>4;
}


template <typename Point, typename Vector>
bool checkSymmetry(const vector<Point> & points, const Vector & axis, const Point& center) {
	typedef Eigen::Matrix<double, Eigen::Dynamic, 4> MatrixXd;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;
	unsigned int size = points.size();
	MatrixXd A(size, 3);
	
	for (int i = 0; i < size; i++) {
		A(i, 0) = (double)points[i][0]*1.0;
		A(i, 1) = (double)points[i][1]*1.0;
		A(i, 2) = (double)points[i][2]*1.0;
		A(i, 3) = 1.0;
	}
	set<Point> initialSet(points.begin(), points.end());
	for (double angle = M_PI/4; angle <= 2 * M_PI; angle+=M_PI/4) {
		set<Point> rotatedSetOfPoints;
		Eigen::Affine3d rot(Eigen::AngleAxisd(angle, Eigen::Vector3d(axis[0], axis[1], axis[2])));
		Eigen::Matrix<double, Eigen::Dynamic, 4> m = A * rot.matrix();
		for (int i = 0; i < size; i++) {
			if (m(i, 0) != 0 || m(i, 1) != 0 || m(i, 2) != 0) {
				int posx = m(i, 0);
				int posy = m(i, 1);
				int posz = m(i, 2);
				Point p{posx, posy, posz};
				rotatedSetOfPoints.insert(p);
			}
		}
		if (areAlmostSimilar(initialSet, rotatedSetOfPoints)) {
			return true;
		}
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////
//////////////// MAIN ///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (skeleton)")
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
	int thresholdMin = vm["thresholdMin"].as<int>();
	int thresholdMax = vm["thresholdMax"].as<int>();
	
	
	
	typedef Z3i::Space Space;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef functors::IntervalForegroundPredicate<Image> ThresholdedImage;
	typedef Z3i::KSpace KSpace;
	typedef DGtal::PointVector<3, double> Vector3D;
	typedef DGtal::PointVector<3, int> Point3D;
	typedef LightImplicitDigitalSurface<KSpace, Z3i::DigitalSet >   MyDigitalSurfaceContainer;
	typedef SurfacePoint<Z3i::Point, SCell> SurfacePoint;
	typedef DigitalSurface<MyDigitalSurfaceContainer> MyDigitalSurface;
	typedef WeightedDigitalSurface<MyDigitalSurfaceContainer, SurfacePoint> MyWeightedDigitalSurface;
	typedef MetricAdjacency<Space, 3> Graph;
	typedef BreadthFirstVisitor<MyDigitalSurface> MyBreadthFirstVisitor;
	typedef BreadthFirstVisitor<MyWeightedDigitalSurface> MyWeightedBreadthFirstVisitor;
	typedef BreadthFirstVisitor<Graph, std::set<Z3i::Point> > PointBreadthFirstVisitor;
	typedef MyBreadthFirstVisitor::Node MyNode;
	typedef MyBreadthFirstVisitor::Size MySize;

	typedef RosenProffittLengthEstimator<set<SCell>> LengthEstimator;
	typedef MSTTangent<Point3D> Tangent;
	typedef Pencil<Point3D, Tangent, Vector3D> Pencil;

	
	
	
	trace.beginBlock("Reading file...");
	Image image = VolReader<Image>::importVol(inputFilename);
	trace.endBlock();
	KSpace ks;
	ks.init( image.domain().lowerBound(), 
			 image.domain().upperBound(), true );
	KSpace::SCellSet boundary;
	Z3i::DigitalSet set3d (image.domain());
	SetFromImage<Z3i::DigitalSet>::append<Image> (set3d, image, 
												  thresholdMin, thresholdMax);
	trace.info() << "Finding a bel" << endl;
	Z3i::SCell bel = Surfaces<KSpace>::findABel( ks, set3d, 100000 );
	//bel = Z3i::SCell({465,516,289},false);
	trace.info() << bel << endl;
	typedef SurfelAdjacency<KSpace::dimension> MySurfelAdjacency;
	MySurfelAdjacency surfAdj( true );
	MyDigitalSurfaceContainer* ptrSurfContainer = 
		new MyDigitalSurfaceContainer( ks, set3d, surfAdj, bel );
	MyDigitalSurface digSurf( ptrSurfContainer ); // acquired
	//! [volBreadthFirstTraversal-SetUpDigitalSurface]

    set<SCell> surfelSet;
	map<SCell, Point> surfaceVoxelSet;
	DigitalSet setPredicate(image.domain());
	
	for (auto it = digSurf.begin(), itE = digSurf.end(); it!=itE; ++it) {
		surfelSet.insert(*it);
		setPredicate.insert(ks.sCoords(*it));
	}
	vector<SurfacePoint> weightedSurfaceV = computeSurfelWeight<SurfacePoint>(ks, set3d, surfelSet, surfaceVoxelSet);
	MyWeightedDigitalSurface weightedSurface(digSurf, weightedSurfaceV);
	MyNode node;
	bool isPathFound = false;
	
    trace.beginBlock("Compute path");
	Z3i::SCell end;
	int numberToFind = 1;
	int i = 0;
	vector<Pencil> pencils;
	set<SCell> thePath;


	QApplication application(argc,argv);
	Viewer3D<> viewer( ks);
	viewer.show();
	viewer << CustomColors3D(Color::Green, Color::Green) << bel;
	while (i < numberToFind) {
		bel = Surfaces<KSpace>::findABel( ks, set3d, 100000 );
		bel = {{27,10,115}, false};
		viewer << CustomColors3D(Color::Green, Color::Green) << bel;
		MyBreadthFirstVisitor visitor( digSurf, bel );
		isPathFound = false;
		i++;
		map<Z3i::SCell, Z3i::SCell> aMapPrevious;
		trace.info() << bel << endl;
		int stepForShortestPath = -2;
		double distanceForShortestPath = numeric_limits<double>::max();
		while (!visitor.finished() && !isPathFound) {
			if (node.second > 100) isPathFound = true;
			node = visitor.current();
			vector<Z3i::SCell> neighbors;
			back_insert_iterator<vector<Z3i::SCell> > iter(neighbors);
			visitor.graph().writeNeighbors(iter, node.first);

			for (auto it = neighbors.begin(), itE = neighbors.end(); it != itE; ++it) {
				auto itMapExisting = aMapPrevious.find(*it);
				if (itMapExisting == aMapPrevious.end())
				{
					aMapPrevious[*it] = node.first;
				}
				else {
					Z3i::SCell tmp = node.first;
					Z3i::SCell tmp2 = aMapPrevious[*it];
					vector<Z3i::Point> aDirectionV;
					vector<Z3i::Point> anotherDirectionV;
					set<Z3i::SCell> path1;
				    set<Z3i::SCell> path2;
				
					while (tmp != bel) {
						path1.insert(tmp);
						Z3i::SCell previous = aMapPrevious[tmp];
						aDirectionV.push_back(ks.sCoords(tmp) - ks.sCoords(previous));
						tmp = previous;
					}
			
					while (tmp2 != bel) {
						path2.insert(tmp2);
						Z3i::SCell previous = aMapPrevious[tmp2];
						anotherDirectionV.push_back(ks.sCoords(tmp2) - ks.sCoords(previous));
						tmp2 = previous;
					}
				    set<Z3i::SCell> x;
					set_union(path1.begin(), path1.end(), path2.begin(), path2.end(), inserter(x, x.end()));
							  
					Z3i::Point center = (ks.sCoords(*it) + ks.sCoords(bel)) / 2;
					//check both paths going back			
					vector<Z3i::SCell> intersection;
					for (auto it = path1.begin(), itE = path1.end(); it!=itE; ++it) {
						for (auto it2 = path2.begin(), itE2 = path2.end(); it2!=itE2; ++it2) {
							if (*it == *it2) {
								intersection.push_back(*it);
							}
						}
					}
					if ((int)node.second == stepForShortestPath + 1) {
						isPathFound = true;
					}
					else if (intersection.size() == 0) {
						vector<Z3i::Point> correspondingPointInPath;
						for (auto it = x.begin(), itE = x.end(); it != itE; ++it) {
							correspondingPointInPath.push_back(surfaceVoxelSet[*it]);
						}
						Vector3D normal = computePlaneNormal<Vector3D>(correspondingPointInPath);
						bool symmetry = checkSymmetry(correspondingPointInPath, normal, center);
						typedef StandardDSS6Computer<vector<Z3i::Point>::iterator,int,8> SegmentComputer;  
						typedef GreedySegmentation<SegmentComputer> Segmentation;
						SegmentComputer algo;
						Segmentation s(correspondingPointInPath.begin(), correspondingPointInPath.end(), algo);
						int cpt = 0; //Cpt corresponds to the number of straight lines in path
						double length = 0;
						for (auto er = s.begin(), erE = s.end(); er!=erE; ++er) {
							SegmentComputer current(*er);
							Z3i::Point direction;
							PointVector<3, double> intercept, thickness;
							current.getParameters(direction, intercept, thickness);
							if (direction != Z3i::Point::zero) {
								cpt++;
								length += euclideanDistance(*(current.end()-1), *(current.begin()));
							}
						}
						if ((int)node.second == stepForShortestPath) {
							viewer << CustomColors3D(Color::Yellow, Color::Yellow) << *it;
							double distance = euclideanDistance(surfaceVoxelSet[*it], surfaceVoxelSet[bel]);
							if(distance < distanceForShortestPath) {
								distanceForShortestPath = distance;
							    thePath = x;
							}
						}
						
						else if (symmetry && stepForShortestPath == -2) {
							stepForShortestPath = node.second;
							viewer << CustomColors3D(Color::Red, Color::Red) << center;
							end = *it;
							//Visualize both paths
							viewer << CustomColors3D(Color::Yellow, Color::Yellow) << *it;
							LengthEstimator l;
							//double distance = l.eval(x.begin(), x.end());
							double distance = euclideanDistance(surfaceVoxelSet[*it], surfaceVoxelSet[bel]);
							if (distance < distanceForShortestPath) {
								distanceForShortestPath = distance;
								thePath = x;
							}

							pencils.push_back({center, normal});
						}

					}
				}
			}
//			viewer << CustomColors3D(Color::Green, Color::Green) << node.first;
			visitor.expand();
		}
		for (auto it = thePath.begin(), itE = thePath.end(); it != itE; ++it) {
			viewer << CustomColors3D(Color::Red, Color::Red) << *it;
		}
		auto visitedV = visitor.visitedVertices();
		for (auto it = visitedV.begin(), itE = visitedV.end();
			 it != itE; ++it) {
			viewer << CustomColors3D(Color::White, Color::White) << *it;
		}
	}
	trace.endBlock();

	
//	SliceUtils::slicesFromPlanes(viewer, pencils, image, "img/slice");
	const Color  CURVE3D_COLOR( 100, 100, 140, 128 );
	for (auto it = set3d.begin(); it != set3d.end(); ++it) {
		viewer <<CustomColors3D(CURVE3D_COLOR, CURVE3D_COLOR)<< (*it);
	}  
	
	viewer << Viewer3D<>::updateDisplay;
	application.exec();

	return 0;
}
