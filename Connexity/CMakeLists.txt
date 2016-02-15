CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(Connexity)
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
endif( Boost_FOUND )

include_directories(../include)
ADD_EXECUTABLE(MakeConnexGeometry MakeConnexGeometry.cpp)
TARGET_LINK_LIBRARIES(MakeConnexGeometry ${DGTAL_LIBRARIES} ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})
ADD_EXECUTABLE(MakeConnex MakeConnex.cpp)
TARGET_LINK_LIBRARIES(MakeConnex ${DGTAL_LIBRARIES} ${Boost_LIBRAIRIES} ${Boost_PROGRAM_OPTIONS_LIBRARY})