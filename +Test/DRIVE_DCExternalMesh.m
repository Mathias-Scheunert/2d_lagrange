% Script for setting up the 2.5D DC problem.
%
% Problem in 3D:
%      x = [x, y, z]
%   \phi = \phi(x)
%
%   -\div(\sigma\grad(\phi)) = I \dirac(x_0) in Omega
%                       \phi = 0             at d_Omega_1 (earth interior)
%               d_\phi / d_n = 0             at d_Omega_2 (surface)
%
% Problem in 2.5D (at z = 0!):
%       x = [x, y]
%   \phi' = \phi'(x, k_Z)
%
%   -\div(\sigma\grad(\phi')) + k_z^2 \sigma \phi' = I/2 \dirac(x_0)
%                    \phi' = 0             at d_Omega_1 (earth interior)
%            d_\phi' / d_n = 0             at d_Omega_2 (surface)
%
% such that
%    \phi = 2/pi \int_{0}^{\inf} \phi' dk_z
%
% 2.5D Variational problem:
%   a(u,v,k_z) = \int_Omega \grad(\phi') * \sigma \grad(v) + ...
%                k_z^2 \int_Omega \phi' * \sigma v
%   f(v)        = I \int_{Omega} \dirac(x_0) v
%
% Numerical integration over wavenumbers:
%   \phi = \sum_{l = 1}^{N} w_l \phi'(k_{z,l})
%
% Coordinate system (Kartesian):
%  0,0 ------>
%      |      x
%      |
%      |
%      v
%       y
%
% REMARKS
%   The 2.5D approach uses the bahavior of a known analytic solution
%   (e.g. point source in 3D half-space) to set up the numerical
%   integration of the seperate 2D solutions (including choice of
%   wavenumbers and weights) which forms the 2.5D solution.
%   This can only act as an approximation for the solution of an arbitrary
%   shaped underground!
%   -> I.e. the solutions should be veryfied by comparing them to a full 3D
%   simulation.
%   Furthermore, this implies that only one single point source can be
%   treated by this approach such that multi-pole arrangements have to be
%   simulated by adding up solutions for several sources.
%
%   2.5D approach derived from:    Dey A., Morrison H.F.; 1979
%   num. integration derived from: Ralph-Uwe Börner (pers. communication)
%                                  Bing Z.; 1998 (Dissertation)
%                                  Xu S.; 2000

%% Set up script.

% Clean up and set verbosity.
clean();
warning('on');
verbosity = pick(2, false, true);

%% Set up disctrete DC fwd problem.

% Define type of numerical integration approach.
FT_type = pick(1, 'Boerner', 'Bing', 'Xu');

% Define source and receiver locations at earth's surface.
% (With an arbitrary topography)
TX.type = 'point_exact';
TX.coo = [1, 1];
TX.val = 1;
RX.coo = [0, 0; 0, 5];

% Define mesh.
mesh_type = 'gmsh_load';
file_name = [pwd, '/+Test/', '2D_circ_complex.msh'];

% Set up domain parameter.
params.val = [1e0, 1e-2];
params.name = {'inner_part', 'outer_part'};

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
bnd.val = {{[];  0; []}, ... % 1 for Dirichlet
           { 0; []; 0}}; ... % 2 for Neumann
bnd.name = {'surface'; 'subsurface'; 'tunnel'};
bnd.quad_ord = 1;

%% Set up FEM.

% Define number of grid refinements.
refinement = 1;

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;
fwd_params.param = params;
clear('TX', 'RX', 'bnd', 'FT_type', 'FE_order', ...
      'refinement', 'param');

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, ...
    'name', file_name, ...
    'ref', fwd_params.ref, ...
    'verbosity', verbosity, ...
    'TX', fwd_params.TX.coo, ...
    'RX', fwd_params.RX.coo);

%% Set up parameter vector.

param = Param.initParam(mesh, fwd_params.param);

%% Assemble 2.5D DC problem.

[fe, sol, FT_info] = App_DC.assembleDC25D(mesh, param, fwd_params, verbosity);

%% Solve 2.5D DC problem.

u_FE = App_DC.solveDC25D(fe, sol, FT_info, verbosity);

%% Visualize solution.

Plot.plotMesh(mesh, param);
colorbar;
Plot.plotSolution(fe, mesh, sum(u_FE, 2));
Plot.plotSolution(fe, mesh, sum(u_FE, 2), 'style', '3D');
