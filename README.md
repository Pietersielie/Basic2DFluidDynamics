![An animation of a tesla valve in the simulated wind tunnel](https://github.com/Pietersielie/Basic2DFluidDynamics/blob/main/samples/TVAColour.gif)
# Basic2DFluidDynamics
A basic wind tunnel simulation, modelling user provided 2D shapes for a given speed. Only supported on Unix, though nothing prevents it from working on Windows as well with appropriate compilers. This implementation is based on the implementation described by Robert Bridson in the [Fluid simulation: SIGGRAPH 2007 course notes](https://dl.acm.org/doi/pdf/10.1145/1281500.1281681).

## Dependencies
1. [CMake](https://cmake.org/) (used for compiling and dealing with the project dependencies).
2. [Eigen C++ library](https://eigen.tuxfamily.org/index.php?title=Main_Page) (used for solving the Poisson equation with the Eigen conjugate gradient solver). Make sure to change line 8 in [CMakeLists.txt](https://github.com/Pietersielie/Basic2DFluidDynamics/blob/main/CMakeLists.txt) to point to your Eigen installation.
3. [OpenCV](https://opencv.org/) (used for reading the mask images and writing the output images).

Eigen and OpenCV can be installed in whatever way you wish, as long as CMake can find them.

## Compilation process
```
git clone https://github.com/Pietersielie/Basic2DFluidDynamics.git
cd Basic2DFluidDynamics
mkdir build && mkdir ImageOutput
cd build
cmake ..
make -j4
```

## Usage instructions
Run the executable created (Basic2DFluidDynamics/build/Basic2DFluidDynamics).
The program has 4 optional arguments: `./Basic2DFluidDynamics [shapeFileName] [WindSpeed] [simulationTime] [interval]`

- **[shapeFileName]**: The path to a file containing the shape you want to model (see samples in the samples folder) 
- **[WindSpeed]**: The speed that the wind tunnel should model. If greater than 0, the wind is from left to right, if negative, it is from right to left.
- **[simulationTime]**: How many seconds the simulation should run for (in simulation time, actual time depends on the hardware and size of the wind tunnel).
- **[interval]**: How many time steps between each image being output to the ImageOutput folder. Note that this is simulation ticks, not actual time.

Take note that the cell size is set to 20 metres on a side -- change this in line 35 of `Basic2DFluidDynamics.cpp` if you desire. The simulation timestep is automatically estimated to be the maximum that still conforms to the $\Delta t \leq \frac{5\Delta x}{u_{max}}$ stability limit recommended by Bridson.
