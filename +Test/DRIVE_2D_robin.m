% Script for solving 2D Laplace problem for manufactured solution.
%
% Coordinate system (Kartesian):
%  0,0 ------>
%      |      x
%      |
%      |
%      v
%       y
%
% COPYRIGHT
%   The script uses the suggested approach from J. Blechta (curlcurl repo).

%% Set up script.

% Clean up and set verbosity.
clean();
warning('on');
verbosity = pick(1, false, true);

% Define number of uniform grid refinements.
refinement = 2;

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

%% Set up disctrete DC fwd problem.

% Define background conductivity.
param_info = struct();
param_info.val = 1;
param_info.name = {'entire'};

% Define mesh.
% Note: TX/RX positions may not be part of the vertex list of 'cube' mesh.
mesh_type = 'disc';

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dtn', 'neumann'};
%         ymin       ymax       xmin       xmax
% @(x, n) 1/pi
dtn_fun = struct();
dtn_fun.f = @(k) @(x, y, n) 1/pi;
dtn_fun.quad_ord = 2 * FE_order;
bnd.val = {{[]; dtn_fun.f}, ...  % 1 for D-t-N operator
           { 0;        []}};     % 2 for Neumann
bnd.name = {'surface', 'subsurface'};
bnd.quad_ord = dtn_fun.quad_ord;

% Summarize parameter.
fwd_params = struct();
fwd_params.TX.coo = [0, 0];
fwd_params.TX.val = 1;
fwd_params.TX.type = 'point_exact';
fwd_params.RX.coo = [];
fwd_params.bnd = bnd;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, ...
                'ref', fwd_params.ref, 'verbosity', verbosity, ...
                'TX', [0, 0], ...
                'dom_name', param_info.name{:});

%% Set up parameter vector.

% Set up parameter vector.
param = Param.initParam(mesh, param_info);

%% Assemble 2D DC problem.

[fe, sol] = App_DC.assembleDC2D(mesh, param, fwd_params, verbosity);

%% Solve 2D DC problem.

% Numeric.
u_FE = FeL.solveFwd(sol, fe, verbosity);

% Manufactured.
ref_fun = struct();
ref_fun.f = @(x, y) 1-1/pi * log(norm([x; y]));
ref_fun.quad_ord = dtn_fun.quad_ord + 1;
u_ref = arrayfun(@(x, y) ref_fun.f(x, y), fe.DOF_maps.DOF_coo(:, 1), ...
                                          fe.DOF_maps.DOF_coo(:, 2));
% Compute symbolic grad of exact solution.
s_x = sym('x', 'real');
s_y = sym('y', 'real');
u_ref_grad_sym = gradient(ref_fun.f(s_x, s_y), [s_x, s_y]);
ref_fun.grad = matlabFunction(u_ref_grad_sym, 'Vars', {s_x, s_y});

% Obtain potentials at RX positions.
phi_FE = fe.I * u_FE;

% Calculate error.
[err_L2, err_H1] = Test.getError(mesh, fe, u_FE, ref_fun);
fprintf('Err L2: %f, Err H1: %f \n', err_L2, err_H1);

%% Visualize.

fig1 = figure(1);
set(fig1, 'Units', 'Normalized', 'OuterPosition', [0.25, 0, 0.5, 1]);
subplot(2, 1, 1);
Plot.plotSolution(fe, mesh, u_ref, 'style', '3D');
title('reference');
zl = zlim;
subplot(2, 1, 2);
Plot.plotSolution(fe, mesh, u_FE, 'style', '3D');
title('FE');
zlim(zl);
