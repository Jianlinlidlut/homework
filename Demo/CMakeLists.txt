cmake_minimum_required(VERSION 3.6)

# Name the project
project(MyDemo)

# Find MeshLib
find_package(MeshLib REQUIRED)
if (NOT MeshLib_FOUND)
    message(FATAL_ERROR "MeshLib Not Found!")
endif (NOT MeshLib_FOUND)

# Find freeglut
find_package(freeglut REQUIRED)
if (NOT freeglut_FOUND)
    message(FATAL_ERROR "freeglut Not Found!")
endif (NOT freeglut_FOUND)

# Include MeshLib core and the 'include' of the project itsself
include_directories(
    "include"
    "${MeshLib_DIR}/core"
    "${freeglut_DIR}/include")

# Manually add the sources
set(SOURCES_ONE ${MeshLib_DIR}/core/bmp/RgbImage.cpp)

# Also, the file(GLOB...) allows for wildcard additions
file(GLOB SOURCES_TWO 
    "include/*.h"
    "src/*.cpp")

# Add an executable target called MyDemo to be build from 
# the source files.
add_executable(MyDemo ${SOURCES_ONE} ${SOURCES_TWO})

# Link the libraries of freeglut
target_link_libraries(MyDemo 
    "${freeglut_DIR}/lib/x64/freeglut.lib")

# Install the executeable program in the bin folder of the
# current project directory.
install(TARGETS MyDemo DESTINATION ${PROJECT_SOURCE_DIR}/bin)
