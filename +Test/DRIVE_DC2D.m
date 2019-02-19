% Script for setting up the 2D DC problem.
%
% Considered is a homogenous electrical background field in a rectangular
% domain including a infinit long conductive cylinder, perpendicular to 
% the E-field.
%
% Problem in 2D.
%      x = [x, y]
%   \phi = \phi(x)
%           
%   -\div(\sigma\grad(\phi)) = I \dirac(x_0) in Omega
%                       \phi = phi_1         at d_Omega_1 (left)
%                       \phi = phi_2         at d_Omega_2 (right)
%               d_\phi / d_n = 0             at d_Omega_3 (top, bottom)
%
% 2D Variational problem:
%   a(u,v) = \int_Omega \grad(\phi') * \sigma \grad(v) + ...
%                \int_Omega \phi' * \sigma v
%   f(v)   = 0
%
% Coordinate system (Kartesian):
%  0,0 ------> 
%      |      x
%      |
%      |
%      v
%       y

%% Set up script.

% Clean up and set verbosity.
clean();
warning('on');
verbosity = pick(2, false, true);

%% Set up mesh.

mesh = Mesh.initMesh('gmsh_load', 'name', '+Test/2D_cylinder.msh');

%% Set up parameter vector.

info = struct();
info.name =  {'background', 'anomaly'};
% For testing.
% info.val = 1./[1e0,          1e0];
% info.val = 1./[1e0,          1e3];
info.val = 1./[1e3,          1e-3];
param = Param.initParam(mesh, info);

%% Set up disctrete DC fwd problem.

% Define source and receiver locations at earth's surface.
% (With an arbitrary topography)
RX.coo = [];
% For testing (only for two opposing sources, hom. N-BC are meaningful).
% TX.type = 'point_exact';
% TX.val = [1;-1];
% TX.coo = [-1.75, 0; ...
%            1.75, 0];
%
TX.type = 'none';

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
% For testing.
% bnd.val = {{[];    [];   [];  []}, ...  % 1 for Dirichlet
%            {0;     0;    0;    0}}; ... % 2 for Neumann
%
bnd.val = {{[];      [];   10;     -10}, ... % 1 for Dirichlet
           {0;     0;      [];     []}}; ... % 2 for Neumann
bnd.name = {'top'; 'bot'; 'left'; 'right'};
bnd.quad_ord = 1;

%% Set up FEM.

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FE_order = FE_order;
clear('TX', 'RX', 'bnd', 'FE_order');
  
%% Assemble 2D DC problem.

[fe, sol] = App_DC.assembleDC2D(mesh, param, fwd_params, verbosity);

%% Solve 2D DC problem.

u_FE = FeL.solveFwd(sol, fe, verbosity);

%% Visualize solution.

Plot.plotGradient(fe, mesh, u_FE(:,1), 'sign', 'neg', 'param', param);
title('Normalized E-field for cylinder in a homogeneous media.');

% Plot.plotGradient(fe, mesh, u_FE(:,2), 'sign', 'neg', 'param', param);
% title('Normalized E-field for cylinder in a homogeneous media.');

% Note: in order to obtain a multipol solution, the solution vectors are
% just summed up.
% Plot.plotGradient(fe, mesh, sum(u_FE, 2), 'sign', 'neg', 'param', param);
% title('Normalized E-field for cylinder in a homogeneous media.');

Plot.plotSolution(fe, mesh, sum(u_FE, 2));
title('Electrical potential for cylinder in a homogeneous media.');