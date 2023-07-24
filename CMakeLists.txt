cmake_minimum_required(VERSION 3.5)

project(VolcanoWeather LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -march=native -fopenmp -fPIC -I /usr/include/eigen3/")

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(lib/)
include_directories(src)

FILE(GLOB_RECURSE HEADER_FILES lib/*.h lib/*.hpp src/*.h src/*.hpp)

add_executable(FlatMACTest
    src/abstractEngine.h
    src/advectionEngine.h
    src/obstaclesEngine.h
    src/pressureEngine.h

    src/airLayer.h
    src/planet.h

    src/AbstractPWMDataStructure.h
    src/Coordinate.h
    src/flatStaggeredGrid.h

    src/atmofuncs.h
    src/mathUtils.h
    src/vizUtils.h

    src/FlatMACTest.cpp
    src/abstractEngine.cpp
    src/planet.cpp
    src/Coordinate.cpp
)