#include <iostream>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/writers/VolWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/images/ImageSelector.h"
#include "DGtal/io/writers/ITKWriter.h"

using namespace std;
using namespace DGtal;
namespace po = boost::program_options;

int main(int argc, char** argv) {

	po::options_description general_opt("Allowed options are: ");
	general_opt.add_options()
		("help,h", "display this message")
		("input,i", po::value<std::string>(), "vol file (.vol) , pgm3d (.p3d or .pgm3d, pgm (with 3 dims)) file or sdp (sequence of discrete points)" );

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
	typedef Z3i::Space Space;
	typedef Z3i::KSpace KSpace;
	typedef HyperRectDomain<Space> Domain;
	typedef ImageSelector<Domain, unsigned char>::Type Image;
	Image image = GenericReader<Image>::import(inputFilename);
	Domain domain = image.domain();
	Z3i::Point lower = image.domain().upperBound();
	int lx = lower[0], ly = lower[1], lz = lower[2];

	size_t lastindex = inputFilename.find_last_of(".");
	string rawname = inputFilename.substr(0, lastindex);
	string outputFilename = rawname + "." + to_string(lx) + "x" + to_string(ly) + "x" + to_string(lz) + ".vol";

	std::ofstream out;

	out.open(outputFilename);

	//We scan the domain instead of the image because we cannot
	//trust the image container Iterator
	for(typename Image::Domain::ConstIterator it = domain.begin(), itend=domain.end();
	    it!=itend;
	    ++it)
	{
		Image::Value val = image( (*it) );
		out.put(val);
	}

	out.close();



	return 0;
}
