% Script for setting up the 2D Poisson problem.
%
% Problem:
%   -nabla²(u) = [0,f] in Omega = ([-1,1] x [-1,1])
%       u(x,y) = g     at d_Omega
%
% Variational problem:
%   a(u,v) = int_Omega grad_u grad_v d(x,y)
%   f(v)   = int_Omega f v           d(x,y)
%
% Galerkin discretization:
%   sum_{i = 1:num(DOF)} u_i int_Omega grad_u_i grad_v_j d(x,y) 
% forall
%   j = 1:num(DOF)

%% Set up script.

% If required, update path using startup.m from the submodule toolbox.
path_req = {'generic'};
path_req = cellfun(@(x) {[pwd, '/', x]}, path_req);
path_cur = regexp(path, regexptranslate('escape', ':'), 'split');
path_missing = ~any(cell2mat(cellfun(@(x) strcmp(x, path_req), path_cur, ...
        'UniformOutput', false).'));
if any(path_missing)
    addpath(path_req{:});
end

% Clean up and set verbosity.
kill();
warning('on');
debug = pick(1, false, true);
verbosity = pick(2, false, true);

%% Set up disctrete Laplace fwd problem.

if debug
    profile on
end
if verbosity
   fprintf('2D Laplace/Poisson problem\n'); 
end

% Define (geometry) scaling factor.
scale = pick(1, 1, 50);

% Define observation points.
RX = pick(1, ...
    [linspace(-0.9, 0.9, scale * 15).', ...
    linspace(0.9, -0.9, scale * 15).'], ... % diagonal profile
    [zeros(scale * 15, 1), ...
    linspace(-0.9, 0.9, scale * 15).']);    % axis parallel profile
RX = scale * RX;

% Define source point and strength.
[Tx, Txp, Txh] = deal(struct());
TXp.type = 'point';
TXp.coo = scale * pick(2, [.1, .3], [0, 0]);
TXp.val = 1; % discrete    Poisson problem
TXp.ref_sol = RefSol.getPeakFunction(TXp.val, TXp.coo);
TXh.type = 'reference';
TXh.val = 1; % homogeneous Poisson problem
TXh.ref_sol = RefSol.getConstFunction(TXh.val);
%
TX = pick(1, TXp, TXh);

% Define boundary conditions.
% ([in]homogeneous Dirichlet conditions,pequal for all four sides)
bnd = struct();
bnd.type = 'dirichlet';
%          bot top left right
bnd.val = [.1,  .1,  .3,   .3    ].';
%
bnd.val = pick(2, ...
    bnd.val * 0, ... %   homogeneous DRB
    bnd.val);        % inhomogeneous DRB

% Define outermost grid boundaries.
x = scale * [-1, 1];
y = scale * [-1, 1];

% Choose basic grid type.
mesh_type = pick(2, 'rhomb', 'cube');

% Set number of grid refinements.
ref_steps = 4;

% Set up order of Lagrange elements.
order = pick(2, 1, 2);

% Print status.
if verbosity
   fprintf(sprintf('- use "%s" basic mesh\n', mesh_type));
   fprintf(sprintf('- use "%d" mesh refinements\n', ref_steps));
   fprintf(sprintf('- use "%s" source\n', TX.type));
   if all(bnd.val) == 0
       fprintf(sprintf('- use homogeneous "%s" boundary conditions\n', bnd.type));
   else
       fprintf(sprintf('- use inhomogeneous "%s" boundary conditions\n', bnd.type));
   end
   fprintf(sprintf('- use oder "%d" Lagrange elements\n', order));
end
if verbosity
   fprintf('... FE-FWP set up.\n \n'); 
end

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, [x, y], ref_steps, verbosity);

%% Set up Parameter.

% Set disturbed area (equals vertical dike).
x_dist = scale * [0, 1];
y_dist = scale * [0.25, 0.35];

% Define parameter of cunductivity (conductor within resistor).
back = 1/100;
dist = pick(1, back, 20);

% Find all cell midpoints using barycentric coordinates.
lambda_mid = 1/3 + zeros(3, 1);
cell_mid = cellfun(@(x) ...
    {[lambda_mid(1)*x(1,1) + lambda_mid(2)*x(2,1) + lambda_mid(3)*x(3,1); ...
    lambda_mid(1)*x(1,2) + lambda_mid(2)*x(2,2) + lambda_mid(3)*x(3,2)]}, ...
    mesh.cell2cord);

% Find cells, belonging to distrubed area.
cell_dist = cellfun(@(x) (x(1) >= x_dist(1) && x(1) <= x_dist(2)) && ...
        (x(2) >= y_dist(1) && x(2) <= y_dist(2)), cell_mid);
    
% Set background parameter for grid.
param = back + zeros(length(mesh.cell2vtx), 1);

% Set parameter for disturbed area.
param(cell_dist) = dist;

%% Set up FE structure.

fe = Fe.initFiniteElement(order, mesh, RX, verbosity);

%% Set up FEM linear System.

% Set up system matrix.
% (for Poisson/Laplace, this only comprises the stiffness matrix)
sol.A = Fe.assembleStiff(fe, param, verbosity);

% Set up rhs vector.
sol.b = Fe.assembleRHS(fe, mesh, TX, verbosity);

% Handle boundary conditions.
[sol, bnd] = Fe.treatDirichlet(fe, mesh, sol, bnd, verbosity);
if verbosity
   fprintf('... Linear system and BC set up.\n \n'); 
end

%% Solve fwd problem.

% Get solution at DOF.
u = Fe.solveFwd(sol, fe, verbosity);

%% Plot solution.

% Solution field.
Plot.plotSolution(fe, mesh, u, param, verbosity);

% Add profile.
% Get solution at RX positions.
phi = fe.I * u;
hold on
    plot3(RX(:,1), RX(:,2), phi, 'r', 'LineWidth', 2);
hold off
if debug
    profile viewer
end