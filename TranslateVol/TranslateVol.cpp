#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/readers/VolReader.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/images/ImageHelper.h"
#include "DGtal/graph/DepthFirstVisitor.h"
#include "DGtal/graph/GraphVisitorRange.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/images/imagesSetsUtils/ImageFromSet.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" )
		("output,o",  po::value<std::string>(), "output itk file" )
		("tx,x", po::value<double>()->default_value(0), "translation x")
		("ty,y", po::value<double>()->default_value(0), "translation y")
		("tz,z", po::value<double>()->default_value(0), "translation z")
		("sx,a", po::value<double>()->default_value(1), "scaling x")
		("sy,b", po::value<double>()->default_value(1), "scaling y")
		("sz,c", po::value<double>()->default_value(1), "scaling z")
		("padding,p", po::value<double>()->default_value(0), "padding");


	bool parseOK=true;
	po::variables_map vm;
	try{
		po::store(po::parse_command_line(argc, argv, general_opt), vm);
	}catch(const std::exception& ex){
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
	string outputFilename = vm["output"].as<std::string>();
	double tx = vm["tx"].as<double>();
	double ty = vm["ty"].as<double>();
	double tz = vm["tz"].as<double>();
	double sx = vm["sx"].as<double>();
	double sy = vm["sy"].as<double>();
	double sz = vm["sz"].as<double>();
	double padding = vm["padding"].as<double>();

	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	typedef Z3i::Object26_6 Object;
	Image image = VolReader<Image>::importVol(inputFilename);
	// Z3i::DigitalSet aSet(image.domain());
	// SetFromImage<Z3i::DigitalSet>::append<Image> (aSet, image,
	// 											  0, 255);
	// trace.info() << aSet.size() << endl;
	// Z3i::DigitalSet cleanSet(image.domain());
	// Object obj(Z3i::dt26_6, aSet);
	// trace.info() << obj.computeConnectedness() << endl;
	// Z3i::DigitalSet & S = obj.pointSet();
	// for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
	// 	if (obj.isSimple(*it)) {
	// 	    S.erase(*it);
	// 	}
	// }

	// for (auto it = S.begin(), ite = S.end(); it != ite; ++it) {
	// 	cleanSet.insert(*it);
	// }
	// Image anotherImage(image.domain());
	// imageFromRangeAndValue(cleanSet.begin(), cleanSet.end(), anotherImage, 1);

	Domain inputDomain = image.domain();
	Z3i::RealVector translationVector(tx, ty, tz);
	Z3i::Point lowerBound = inputDomain.lowerBound() + translationVector;
	Z3i::Point tmp = (inputDomain.upperBound() - inputDomain.lowerBound());
	Z3i::Point upperBound = Z3i::Point(tmp[0] * sx, tmp[1] * sy, tmp[2] * sz) + lowerBound;

	Z3i::Point paddingPoint(padding, padding, padding);
	Z3i::Domain domain(lowerBound-paddingPoint, upperBound+paddingPoint);
	Image out(domain);
	for (auto it = image.domain().begin(), ite = image.domain().end();
		 it!=ite; ++it) {
		Z3i::Point p = *it;
		if (image(p)>0) {
			for (float i = 0; i < sx; i++) {
				for (float j = 0; j < sy; j++) {
					for (float k = 0; k < sz; k++) {
						Z3i::Point tmp = p - inputDomain.lowerBound();
						Z3i::Point newP = Z3i::Point(tmp[0] * sx + i,
											 		 tmp[1] * sy + j,
													 tmp[2] * sz + k)
							+ lowerBound;
						out.setValue(newP, 255);

					}
				}
			}
		}
	}


	VolWriter<Image>::exportVol(outputFilename, out);
	return 0;
}
