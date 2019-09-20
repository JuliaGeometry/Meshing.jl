# Meshing.jl

This package provides a suite of meshing (isosurface extraction) algorithms.

Algorithms included:
* [Marching Tetrahedra](https://en.wikipedia.org/wiki/Marching_tetrahedra)
* [Marching Cubes](https://en.wikipedia.org/wiki/Marching_cubes)
* [Naive Surface Nets](https://0fps.net/2012/07/12/smooth-voxel-terrain-part-2/)

## What is isosurface extraction?

Isosurface extraction is a common technique to visualize and analyze scalar fields and functions.
It takes scalar data that is structured in grids and converts this into a series of points and faces suitable for 3D visualization, GPU computing, 3D Printing, or other analyses.

There are several applications for these techniques such as:
- Medical Imaging of CT and MRI data
- Visualization of Functions
- Solid Modeling
- Terrain Generation
