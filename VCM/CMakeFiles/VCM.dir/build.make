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
CMAKE_SOURCE_DIR = /home/florent/Documents/DGtal/VCM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/florent/Documents/DGtal/VCM

# Include any dependencies generated for this target.
include CMakeFiles/VCM.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/VCM.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/VCM.dir/flags.make

CMakeFiles/VCM.dir/VCM.cpp.o: CMakeFiles/VCM.dir/flags.make
CMakeFiles/VCM.dir/VCM.cpp.o: VCM.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/florent/Documents/DGtal/VCM/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/VCM.dir/VCM.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/VCM.dir/VCM.cpp.o -c /home/florent/Documents/DGtal/VCM/VCM.cpp

CMakeFiles/VCM.dir/VCM.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/VCM.dir/VCM.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/florent/Documents/DGtal/VCM/VCM.cpp > CMakeFiles/VCM.dir/VCM.cpp.i

CMakeFiles/VCM.dir/VCM.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/VCM.dir/VCM.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/florent/Documents/DGtal/VCM/VCM.cpp -o CMakeFiles/VCM.dir/VCM.cpp.s

CMakeFiles/VCM.dir/VCM.cpp.o.requires:
.PHONY : CMakeFiles/VCM.dir/VCM.cpp.o.requires

CMakeFiles/VCM.dir/VCM.cpp.o.provides: CMakeFiles/VCM.dir/VCM.cpp.o.requires
	$(MAKE) -f CMakeFiles/VCM.dir/build.make CMakeFiles/VCM.dir/VCM.cpp.o.provides.build
.PHONY : CMakeFiles/VCM.dir/VCM.cpp.o.provides

CMakeFiles/VCM.dir/VCM.cpp.o.provides.build: CMakeFiles/VCM.dir/VCM.cpp.o

# Object files for target VCM
VCM_OBJECTS = \
"CMakeFiles/VCM.dir/VCM.cpp.o"

# External object files for target VCM
VCM_EXTERNAL_OBJECTS =

VCM: CMakeFiles/VCM.dir/VCM.cpp.o
VCM: CMakeFiles/VCM.dir/build.make
VCM: /usr/local/lib/libDGtal.so
VCM: /usr/local/lib/libDGtalIO.so
VCM: /usr/lib/x86_64-linux-gnu/libboost_program_options.a
VCM: /usr/local/lib/libDGtal.so
VCM: /usr/lib/x86_64-linux-gnu/libcairo.so
VCM: /usr/lib/x86_64-linux-gnu/libQtOpenGL.so
VCM: /usr/lib/x86_64-linux-gnu/libQtGui.so
VCM: /usr/lib/x86_64-linux-gnu/libQtXml.so
VCM: /usr/lib/x86_64-linux-gnu/libQtCore.so
VCM: /usr/lib/x86_64-linux-gnu/libQGLViewer.so
VCM: /usr/lib/x86_64-linux-gnu/libGLU.so
VCM: /usr/lib/x86_64-linux-gnu/libGL.so
VCM: /usr/lib/x86_64-linux-gnu/libSM.so
VCM: /usr/lib/x86_64-linux-gnu/libICE.so
VCM: /usr/lib/x86_64-linux-gnu/libX11.so
VCM: /usr/lib/x86_64-linux-gnu/libXext.so
VCM: CMakeFiles/VCM.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable VCM"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/VCM.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/VCM.dir/build: VCM
.PHONY : CMakeFiles/VCM.dir/build

CMakeFiles/VCM.dir/requires: CMakeFiles/VCM.dir/VCM.cpp.o.requires
.PHONY : CMakeFiles/VCM.dir/requires

CMakeFiles/VCM.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/VCM.dir/cmake_clean.cmake
.PHONY : CMakeFiles/VCM.dir/clean

CMakeFiles/VCM.dir/depend:
	cd /home/florent/Documents/DGtal/VCM && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/florent/Documents/DGtal/VCM /home/florent/Documents/DGtal/VCM /home/florent/Documents/DGtal/VCM /home/florent/Documents/DGtal/VCM /home/florent/Documents/DGtal/VCM/CMakeFiles/VCM.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/VCM.dir/depend

