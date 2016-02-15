CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Curve)
FIND_PACKAGE(DGtal REQUIRED)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


SET(Boost_USE_STATIC_LIBS   ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)
SET(Boost_FOUND FALSE)
FIND_PACKAGE(Boost 1.46.0 REQUIRED COMPONENTS program_options)
if ( Boost_FOUND )
  message(STATUS "Boost and boost_program_options found.")
  include_directories( ${Boost_INCLUDE_DIRS} )
  SET(DGtalToolsLibDependencies ${DGtalToolsLibDependencies} 
     ${Boost_LIBRAIRIES}  
     ${Boost_PROGRAM_OPTIONS_LIBRARY})
   SET(DGtalLibInc ${Boost_INCLUDE_DIRS})
   message(STATUS ${Boost_INCLUDE_DIRS})
   message(STATUS ${Boost_LIBRAIRIES})  
endif( Boost_FOUND )


find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )
include_directories(../include)
ADD_EXECUTABLE(Curve Curve.cpp)
ADD_EXECUTABLE(Noise testNoise.cpp)
TARGET_LINK_LIBRARIES(Curve ${DGTAL_LIBRARIES})
TARGET_LINK_LIBRARIES(Noise ${DGTAL_LIBRARIES})