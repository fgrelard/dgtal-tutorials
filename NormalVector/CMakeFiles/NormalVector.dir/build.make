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
CMAKE_SOURCE_DIR = /home/florent/Documents/DGtal/NormalVector

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/florent/Documents/DGtal/NormalVector

# Include any dependencies generated for this target.
include CMakeFiles/NormalVector.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/NormalVector.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/NormalVector.dir/flags.make

CMakeFiles/NormalVector.dir/NormalVector.cpp.o: CMakeFiles/NormalVector.dir/flags.make
CMakeFiles/NormalVector.dir/NormalVector.cpp.o: NormalVector.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/florent/Documents/DGtal/NormalVector/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/NormalVector.dir/NormalVector.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/NormalVector.dir/NormalVector.cpp.o -c /home/florent/Documents/DGtal/NormalVector/NormalVector.cpp

CMakeFiles/NormalVector.dir/NormalVector.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/NormalVector.dir/NormalVector.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/florent/Documents/DGtal/NormalVector/NormalVector.cpp > CMakeFiles/NormalVector.dir/NormalVector.cpp.i

CMakeFiles/NormalVector.dir/NormalVector.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/NormalVector.dir/NormalVector.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/florent/Documents/DGtal/NormalVector/NormalVector.cpp -o CMakeFiles/NormalVector.dir/NormalVector.cpp.s

CMakeFiles/NormalVector.dir/NormalVector.cpp.o.requires:
.PHONY : CMakeFiles/NormalVector.dir/NormalVector.cpp.o.requires

CMakeFiles/NormalVector.dir/NormalVector.cpp.o.provides: CMakeFiles/NormalVector.dir/NormalVector.cpp.o.requires
	$(MAKE) -f CMakeFiles/NormalVector.dir/build.make CMakeFiles/NormalVector.dir/NormalVector.cpp.o.provides.build
.PHONY : CMakeFiles/NormalVector.dir/NormalVector.cpp.o.provides

CMakeFiles/NormalVector.dir/NormalVector.cpp.o.provides.build: CMakeFiles/NormalVector.dir/NormalVector.cpp.o

# Object files for target NormalVector
NormalVector_OBJECTS = \
"CMakeFiles/NormalVector.dir/NormalVector.cpp.o"

# External object files for target NormalVector
NormalVector_EXTERNAL_OBJECTS =

NormalVector: CMakeFiles/NormalVector.dir/NormalVector.cpp.o
NormalVector: CMakeFiles/NormalVector.dir/build.make
NormalVector: /usr/local/lib/libDGtal.so
NormalVector: /usr/local/lib/libDGtalIO.so
NormalVector: /usr/lib/x86_64-linux-gnu/libboost_program_options.a
NormalVector: /usr/local/lib/libDGtal.so
NormalVector: /usr/lib/x86_64-linux-gnu/libcairo.so
NormalVector: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
NormalVector: /usr/lib/x86_64-linux-gnu/libQtGui.so
NormalVector: /usr/lib/x86_64-linux-gnu/libQtXml.so
NormalVector: /usr/lib/x86_64-linux-gnu/libQtCore.so
NormalVector: /usr/lib/x86_64-linux-gnu/libQGLViewer.so
NormalVector: /usr/lib/x86_64-linux-gnu/libGLU.so
NormalVector: /usr/lib/x86_64-linux-gnu/libGL.so
NormalVector: /usr/lib/x86_64-linux-gnu/libSM.so
NormalVector: /usr/lib/x86_64-linux-gnu/libICE.so
NormalVector: /usr/lib/x86_64-linux-gnu/libX11.so
NormalVector: /usr/lib/x86_64-linux-gnu/libXext.so
NormalVector: CMakeFiles/NormalVector.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable NormalVector"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/NormalVector.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/NormalVector.dir/build: NormalVector
.PHONY : CMakeFiles/NormalVector.dir/build

CMakeFiles/NormalVector.dir/requires: CMakeFiles/NormalVector.dir/NormalVector.cpp.o.requires
.PHONY : CMakeFiles/NormalVector.dir/requires

CMakeFiles/NormalVector.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/NormalVector.dir/cmake_clean.cmake
.PHONY : CMakeFiles/NormalVector.dir/clean

CMakeFiles/NormalVector.dir/depend:
	cd /home/florent/Documents/DGtal/NormalVector && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/florent/Documents/DGtal/NormalVector /home/florent/Documents/DGtal/NormalVector /home/florent/Documents/DGtal/NormalVector /home/florent/Documents/DGtal/NormalVector /home/florent/Documents/DGtal/NormalVector/CMakeFiles/NormalVector.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/NormalVector.dir/depend

