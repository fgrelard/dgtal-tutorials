# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests

# Include any dependencies generated for this target.
include CMakeFiles/TestSurfaceTraversal.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TestSurfaceTraversal.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TestSurfaceTraversal.dir/flags.make

CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o: CMakeFiles/TestSurfaceTraversal.dir/flags.make
CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o: testSurfaceTraversal.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o -c /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests/testSurfaceTraversal.cpp

CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests/testSurfaceTraversal.cpp > CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.i

CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests/testSurfaceTraversal.cpp -o CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.s

CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.requires:
.PHONY : CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.requires

CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.provides: CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.requires
	$(MAKE) -f CMakeFiles/TestSurfaceTraversal.dir/build.make CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.provides.build
.PHONY : CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.provides

CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.provides.build: CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o

# Object files for target TestSurfaceTraversal
TestSurfaceTraversal_OBJECTS = \
"CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o"

# External object files for target TestSurfaceTraversal
TestSurfaceTraversal_EXTERNAL_OBJECTS =

TestSurfaceTraversal: CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o
TestSurfaceTraversal: CMakeFiles/TestSurfaceTraversal.dir/build.make
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libmpfr.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libgmp.so
TestSurfaceTraversal: /usr/local/lib/libCGAL_Core.so
TestSurfaceTraversal: /usr/local/lib/libCGAL.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libboost_thread.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libboost_system.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libpthread.so
TestSurfaceTraversal: /usr/local/lib/libDGtal.so
TestSurfaceTraversal: /usr/local/lib/libDGtalIO.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libboost_program_options.a
TestSurfaceTraversal: /usr/local/lib/libDGtal.so
TestSurfaceTraversal: /usr/lib/libblas.so
TestSurfaceTraversal: /usr/lib/liblapack.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libmpfr.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libgmpxx.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libgmp.so
TestSurfaceTraversal: /usr/local/lib/libITKIOTransformHDF5-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOTransformInsightLegacy-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOTransformMatlab-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOTransformBase-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKDeprecated-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOBMP-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOBioRad-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOGDCM-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkgdcmMSFF-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkopenjpeg-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkgdcmDICT-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkgdcmIOD-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkgdcmDSED-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkgdcmCommon-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOGIPL-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOJPEG-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOLSM-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOTIFF-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitktiff-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkjpeg-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOMeta-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIONIFTI-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIONRRD-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKNrrdIO-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOPNG-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkpng-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOStimulate-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOVTK-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKLabelMap-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKQuadEdgeMesh-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKBiasCorrection-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKPolynomials-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKBioCell-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKDICOMParser-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOSpatialObjects-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOXML-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKFEM-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKMetaIO-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOSiemens-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOGE-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOIPL-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKKLMRegionGrowing-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKVTK-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKWatersheds-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKSpatialObjects-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKMesh-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKPath-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOMesh-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKgiftiio-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKEXPAT-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKniftiio-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKznz-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOCSV-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOHDF5-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkhdf5_cpp-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkhdf5-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkzlib-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOMRC-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKOptimizersv4-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKOptimizers-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKStatistics-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkNetlibSlatec-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKVideoIO-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKIOImageBase-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKVideoCore-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKCommon-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkdouble-conversion-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitksys-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libITKVNLInstantiation-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkvnl_algo-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkv3p_lsqr-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkvnl-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkvcl-4.6.so.1
TestSurfaceTraversal: /usr/local/lib/libitkv3p_netlib-4.6.so.1
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libcairo.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libQtGui.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libQtXml.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libQtCore.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libQGLViewer.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libGLU.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libGL.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libSM.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libICE.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libX11.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libXext.so
TestSurfaceTraversal: /usr/local/lib/libCGAL_Core.so
TestSurfaceTraversal: /usr/local/lib/libCGAL.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libboost_thread.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libboost_system.so
TestSurfaceTraversal: /usr/lib/x86_64-linux-gnu/libpthread.so
TestSurfaceTraversal: CMakeFiles/TestSurfaceTraversal.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable TestSurfaceTraversal"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TestSurfaceTraversal.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TestSurfaceTraversal.dir/build: TestSurfaceTraversal
.PHONY : CMakeFiles/TestSurfaceTraversal.dir/build

CMakeFiles/TestSurfaceTraversal.dir/requires: CMakeFiles/TestSurfaceTraversal.dir/testSurfaceTraversal.cpp.o.requires
.PHONY : CMakeFiles/TestSurfaceTraversal.dir/requires

CMakeFiles/TestSurfaceTraversal.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TestSurfaceTraversal.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TestSurfaceTraversal.dir/clean

CMakeFiles/TestSurfaceTraversal.dir/depend:
	cd /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests /home/florent/Documents/DGtal/tutorials-copy/include/surface/tests/CMakeFiles/TestSurfaceTraversal.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TestSurfaceTraversal.dir/depend

