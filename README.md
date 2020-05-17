# Final Project: Snow Simulation

This repo is the implementation of my choice of the final project of the open course: GAMES 101 -  Introduction to Computer Graphics. The topic is "A Material Point Method for Snow Simulation" (Stomakhin et al., 2013). The simulation is 3D with the help of the following dependencies:

- OpenMP
- Eigen3
- OpenGL
- GLUT

## Build Instruction

```
mkdir build
cd build
cmake ..
make
```

## Default Scene Description

Owning to the limited time of this final project, the interactive part is yet finished, still, you can try to alter the corresponding lines in main.cpp to change the scene bounding box, and the initial snow particle distribution by filling in spheres, rectangles, and any closed triangle mesh.

The default scene (by 10/May/2020) is "two cows bouncing into each other", and there is also "a very thin snow layer on the floor" to show two cows' trajectory.

## TODO List

Again, due to the deadline (by 10/May/2020), there are many things I would like to implement but cannot finish on time.

A General TODO list is provided below. Any questions related to the fulfillment, please contact elliott.yd.cao@outlook.com. If you have any suggestions or improvements, please file an issue or open a PR directly for simple error fixing.

* Implicit calibration
* Collision detection for the sphere, rectangle, and the explicit triangle mesh
* 3D WordArt snow shape
* Export results with compatible results for Houdini for rendering
* Better multi-threads management for GUI and simulational domain

## Related Links
- For more details of the open course GAMES 101, please refer to the [link](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html).
- For the full copy of the paper, please refer to the [link](https://dl.acm.org/doi/10.1145/2461912.2461948). The copyright holders grant permission to study or personal use without commercial purposes FOR FREE.
- For the collection of my implementations of all assignments of GAMES 101, please refer to this [repository](https://github.com/Eydcao/YanCG).