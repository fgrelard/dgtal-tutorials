# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests

# Include any dependencies generated for this target.
include CMakeFiles/TestGeodesicMetric.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TestGeodesicMetric.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestGeodesicMetric.dir/flags.make

CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o: CMakeFiles/TestGeodesicMetric.dir/flags.make
CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o: testGeodesicMetric.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o -c /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests/testGeodesicMetric.cpp

CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests/testGeodesicMetric.cpp > CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.i

CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests/testGeodesicMetric.cpp -o CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.s

CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.requires:
.PHONY : CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.requires

CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.provides: CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.requires
	$(MAKE) -f CMakeFiles/TestGeodesicMetric.dir/build.make CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.provides.build
.PHONY : CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.provides

CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.provides.build: CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o

# Object files for target TestGeodesicMetric
TestGeodesicMetric_OBJECTS = \
"CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o"

# External object files for target TestGeodesicMetric
TestGeodesicMetric_EXTERNAL_OBJECTS =

TestGeodesicMetric: CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o
TestGeodesicMetric: CMakeFiles/TestGeodesicMetric.dir/build.make
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libmpfr.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libgmp.so
TestGeodesicMetric: /usr/local/lib/libCGAL_Core.so
TestGeodesicMetric: /usr/local/lib/libCGAL.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libboost_thread.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libboost_system.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libpthread.so
TestGeodesicMetric: /usr/local/lib/libDGtalIO.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libboost_program_options.a
TestGeodesicMetric: /usr/local/lib/libDGtal.so
TestGeodesicMetric: /usr/lib/libblas.so
TestGeodesicMetric: /usr/lib/liblapack.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libmpfr.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libgmpxx.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libgmp.so
TestGeodesicMetric: /usr/local/lib/libITKDeprecated-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKDICOMParser-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOCSV-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOHDF5-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOLSM-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOMRC-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOMesh-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKgiftiio-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKOptimizersv4-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKReview-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKLabelMap-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKQuadEdgeMesh-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKPolynomials-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKBiasCorrection-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKBioCell-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOGDCM-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkgdcmMSFF-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkgdcmDICT-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkgdcmIOD-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkgdcmDSED-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkgdcmCommon-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOSpatialObjects-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOXML-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKEXPAT-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKFEM-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKOptimizers-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOBMP-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOBioRad-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOGE-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOGIPL-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOJPEG-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOTIFF-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitktiff-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkjpeg-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOMeta-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKMetaIO-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIONIFTI-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKniftiio-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKznz-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIONRRD-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKNrrdIO-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOPNG-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkpng-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOSiemens-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOIPL-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOStimulate-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOTransformHDF5-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkhdf5_cpp-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkhdf5-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkzlib-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOTransformInsightLegacy-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOTransformMatlab-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOTransformBase-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOVTK-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKKLMRegionGrowing-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkopenjpeg-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKVTK-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKWatersheds-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKStatistics-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkNetlibSlatec-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKSpatialObjects-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKMesh-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKTransform-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKPath-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKVideoIO-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKIOImageBase-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKVideoCore-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKCommon-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkdouble-conversion-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitksys-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libITKVNLInstantiation-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkvnl_algo-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkv3p_lsqr-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkvnl-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkvcl-4.10.so.1
TestGeodesicMetric: /usr/local/lib/libitkv3p_netlib-4.10.so.1
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libcairo.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libQtGui.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libQtXml.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libQtCore.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libQGLViewer.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libGLU.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libGL.so
TestGeodesicMetric: /usr/local/lib/libCGAL_Core.so
TestGeodesicMetric: /usr/local/lib/libCGAL.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libboost_thread.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libboost_system.so
TestGeodesicMetric: /usr/lib/x86_64-linux-gnu/libpthread.so
TestGeodesicMetric: CMakeFiles/TestGeodesicMetric.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable TestGeodesicMetric"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestGeodesicMetric.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestGeodesicMetric.dir/build: TestGeodesicMetric
.PHONY : CMakeFiles/TestGeodesicMetric.dir/build

CMakeFiles/TestGeodesicMetric.dir/requires: CMakeFiles/TestGeodesicMetric.dir/testGeodesicMetric.cpp.o.requires
.PHONY : CMakeFiles/TestGeodesicMetric.dir/requires

CMakeFiles/TestGeodesicMetric.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestGeodesicMetric.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestGeodesicMetric.dir/clean

CMakeFiles/TestGeodesicMetric.dir/depend:
	cd /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests /home/florent/Documents/DGtal/tutorials-copy/include/geometry/tests/CMakeFiles/TestGeodesicMetric.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TestGeodesicMetric.dir/depend

