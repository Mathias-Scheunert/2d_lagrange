% Script for setting up the 2D Poisson problem.
%
% Coordinate system:
%  0,0 ------> 
%      |      x
%      |
%      |
%      v
%       y
%
% Problem:
%   -nabla²(u) = 0 in Omega = ([-1,1] x [-1,1])
%       x(x,y) = g at d_Omega
%
% Variational problem:
%   a(u,v) = 0 = int_Omega grad_u grad_v d(x,y)
%
% Galerkin discretization:
%   sum_{i = 1:num(DOF)} u_i int_Omega grad_u_i grad_v_j d(x,y) 
% forall
%   j = 1:num(DOF)

%% Set up script.

% If required, update path using startup.m from the submodule toolbox.
path_req = {'generic'};
path_req = cellfun(@(x) {[pwd, '/', x]}, path_req);
path_cur = strsplit(path, pathsep);
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
TXp.coo = scale * [.1, .3];
TXp.val = 1; % discrete    Poisson problem
TXh.type = 'homogen';
TXh.val = 1; % homogeneous Poisson problem
TX = pick(1, TXp, TXh);

% Define boundary conditions.
% ([in]homogeneous Dirichlet conditions, equal for all four sides)
bnd = struct();
bnd.type = 'dirichlet';
%          bot top left right
bnd.val = [0,  0,  0,   0    ].';
bnd.val = pick(1, ...
    bnd.val, ...  %   homogeneous DRB
    bnd.val + 5); % inhomogeneous DRB

% Define outermost grid boundaries.
x = scale * [-1, 1];
y = scale * [-1, 1];

% Choose basic grid type.
mesh_type = pick(2, 'rhomb', 'cube');

% Set number of grid refinements.
%                   1  2  3
ref_steps = pick(2, 2, 4, 6);

% Set up order of Lagrange elements.
order = pick(2, 1, 2);

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
   fprintf('... FE-FWP set up.\n ------ \n'); 
end

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, [x, y], ref_steps, verbosity);

if verbosity
   fprintf('... Mesh struct initialized.\n ------ \n'); 
end
%% Set up Parameter.
    
% Set background parameter for grid.
param = 1 + zeros(length(mesh.cell2vtx), 1);

%% Set up FE structure.

fe = Fe.initFiniteElement(order, mesh, RX, verbosity);

if verbosity
   fprintf('... FE struct initialized.\n ------ \n'); 
end

%% Set up FEM linear System.

% Set up system matrix.
% (for Poisson/Laplace, this only comprises the stiffness matrix)
if verbosity
   fprintf('Assemble stiffness matrix ... '); 
end
sol.A = Fe.assembleStiff(fe, param);
if verbosity
   fprintf('done.\n'); 
end

% Set up rhs vector.
if verbosity
   fprintf('Assemble rhs ... '); 
end
sol.b = Fe.assembleRHS(fe, mesh, TX);
if verbosity
   fprintf('done.\n'); 
end

% Handle boundary conditions.
if verbosity
   fprintf('Incorporate Dirichlet BC ... '); 
end
sol = Fe.treatDirichlet(fe, mesh, sol, bnd);
if verbosity
   fprintf('done.\n'); 
end
if verbosity
   fprintf('... Linear system and BC set up.\n ------ \n'); 
end

%% Solve fwd problem.

% Get solution at DOF.
u = Fe.solveFwd(sol, fe, verbosity);

if verbosity
   fprintf('... FWP solved.\n ------ \n'); 
end
%% Plot solution.

% Solution field.
Plot.plotSolution(fe, mesh, u, verbosity);

% Add profile.
% Get solution at RX positions.
phi = fe.I * u;
hold on
    plot3(RX(:,1), RX(:,2), phi, 'r', 'LineWidth', 2);
hold off
if verbosity
   fprintf('... Solution visualized.\n'); 
end
if debug
    profile viewer
end