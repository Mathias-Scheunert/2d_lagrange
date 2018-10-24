# Install

- clone the repository: `git clone git@git.tu-freiberg.de:ng/2d-lagrange.git`

Note that the code was developed at MATLAB 2017b.
<br/>
Hence, there might be some issues w.r.t. to missing MATLAB builtin functions at older version.
<br/>
<br/>
Update: MATLAB 2015 and 2016 supported

# Misc

The motivation for this repo was to get familar with the (implementation of the) FE method.
<br/>
Hence, it provides only the basic functionality which can be extended ad libitum.
<br/>
I used the codestyle, introduced at [(codestyle)](https://git.tu-freiberg.de/ng/toolbox/blob/master/template/codeStyleTemplate.m).

Feel free to use this repo for:
- teaching
- developing (geo)physical applications
- just playing around ... 

# Features

Mesh:
- arbitrary bounded domains (exctracted areas at interior allowed)
- initial mesh variants: uniform splitted rectangle and rectangle domain including a rhombus
- uniform mesh refinement
- Gmsh: 
    - workable Gmsh binary
	- routine to create default DC Gmsh .geo input file (including topography, TX/RX positions)
	- routine to load arbitrary mesh from Gmsh .msh output files (2D and 3D supported)
	- using Gmsh physical tags to identify parameter domains (surfaces) and boundary parts (lines) from .msh

Lagrange elements:
- 2D first and second order

Raviart-Thomas elements:
- 2D zeroth order

Assembling of:
- Lagrange:
    - mass and stiffness matrix
    - tensor representation for derivative of these matrices w.r.t. the parameter
    - interpolation operator for observing solution at arbitrary points within the domain
- Raviart-Thomas:
    - mass and divergence matrix
    - regularization matrix representing smoothness conatraints for unstructured grids
- rhs vector for point source(es)
- rhs vector for given reference function

Boundary conditions:
- homogeneous and inhomogeneous Dirichlet
- homogeneous and inhomogeneous Neumann

Visualization:
- mesh with parameter and edge orientations and normals
- Lagrange FE solution within domain and observations
- Raviart-Thomas elements/basis functions (on local, reference simplex and in global mesh)

# Applications:

The following driver files are available for application:
- DRIVE_DC: 2.5D DC problem (modelling of electrical potential)
    - 1 TX (3D-point source), n RX (arbitrary position in mesh)
    - 3 types of spatial inverse FT
    - in-/homogeneous Dirichlet/Neumann boundary conditions
    - mesh generation via Gmsh (incorporating TX/RX positions, topography)

The following driver files are available for testing:
- DRIVE_TestMesh: Visualization of the underlying unstructured grid and derived information.
- DRIVE_Poisson: FE solution of the Poission equation for homogeneous/point source and in-/homogeneous Dirichlet/Neumann boundary conditions.
- DRIVE_FEvsRefSol: 
    - FE reference solution for an arbitrary polynomial function whose analytic derivative act as rhs of the problem.
    - FE reference solution for a 2D-point source with known analytic solution
    - in-/homogeneous Dirichlet/Neumann boundary conditions
- DRIVE_DCExternal: 
    - 2.5D DC problem with a 3D-point source
    - homogeneous Dirichlet BC at boundaries in subsurface, homogeneous Neumann BC at boundary on surface
    - incorporation of an external mesh, created by Gmsh
- testLagrange:
    - class for matlab unittest runner
    - uses reference solutions from DRIVE_FEvsRefSol and checks convergence rates for different source and BC types
    - provides L2 and H1 error