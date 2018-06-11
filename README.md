# Install

- clone the repository: `git clone git@git.tu-freiberg.de:ng/2d-lagrange.git`

Note that the code was developed at MATLAB 2017b.
<br/>
Hence, there might be some issues w.r.t to missing MATLAB builtin functions at older version.

# Usage

The motivation for this repo was to get familar with the (implementation of the) FE method.
<br/>
Hence, it provides only the basic functionality which can be extended ad libitum.
<br/>
I used the codestyle, introduced at [(codestyle)](https://git.tu-freiberg.de/ng/toolbox/blob/master/template/codeStyleTemplate.m).

There are two driver files available:
- DRIVE_testMesh: Visualization of the underlying unstructured grid and derived information.
- DRIVE_Poisson: FE solution of the Poission equation for homogeneous/point source and in-/homogeneous Dirichlet boundary conditions.

Feel free to use this repo for:
- teaching
- developing (geo)physical appications
- just playing around ...

