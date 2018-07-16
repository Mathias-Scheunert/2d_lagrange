# Install

- clone the repository: `git clone git@git.tu-freiberg.de:ng/2d-lagrange.git`

Note that the code was developed at MATLAB 2017b.
<br/>
Hence, there might be some issues w.r.t. to missing MATLAB builtin functions at older version.
<br/>
<br/>
Update: MATLAB 2015 and 2016 supported

# Usage

The motivation for this repo was to get familar with the (implementation of the) FE method.
<br/>
Hence, it provides only the basic functionality which can be extended ad libitum.
<br/>
I used the codestyle, introduced at [(codestyle)](https://git.tu-freiberg.de/ng/toolbox/blob/master/template/codeStyleTemplate.m).

The following driver files are available:
- DRIVE_TestMesh: Visualization of the underlying unstructured grid and derived information.
- DRIVE_TestCode: 
    - FE solution for an arbitrary polynomial function whose analytic derivative act as rhs of the problem.
    - FE solution for a point source
    - L2 and H1 error convergence test
- DRIVE_Poisson: FE solution of the Poission equation for homogeneous/point source and in-/homogeneous Dirichlet boundary conditions.

Feel free to use this repo for:
- teaching
- developing (geo)physical applications
- just playing around ... 

# Features

Mesh:
- domains bounded by a rectangle
- initial mesh variants: uniform splitted rectangle and arbitrary mesh including a rhombus
- uniform mesh refinement

Lagrange elements:
- 2D first and second order

Assembling of:
- mass and stiffness matrix
- rhs vector for point source(es)
- rhs vector for given reference function
- interpolation operator for observing solution at arbitrary points within the domain

Boundary conditions:
- homogeneous and inhomogeneous Dirichlet
- homogeneous and inhomogeneous Neumann

Visualization:
- mesh with parameter and edge orientations and normals
- FE solution within domain and observations