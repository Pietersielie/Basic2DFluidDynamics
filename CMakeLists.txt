cmake_minimum_required(VERSION 3.5)

project(Basic2DFluidDynamics LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native -fopenmp -fPIC -I /usr/include/eigen3/")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(lib/)
include_directories(src)

FILE(GLOB_RECURSE HEADER_FILES lib/*.h lib/*.hpp src/*.h src/*.hpp)

add_executable(Basic2DFluidDynamics
    src/abstractEngine.h
    src/advectionEngine.h
    src/obstaclesEngine.h
    src/pressureEngine.h

    src/fluidLayer.h

    src/AbstractDataStructure.h
    src/Coordinate.h
    src/flatStaggeredGrid.h

    src/atmofuncs.h
    src/mathUtils.h
    src/vizUtils.h

    src/FlatMACTest.cpp
    src/abstractEngine.cpp
    src/Coordinate.cpp
)
