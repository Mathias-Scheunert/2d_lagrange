% Script for setting up the 2D Poisson problem.
%
% Problem:
%   -\div(\sigma\grad(u)) = ...
%   -\sigma\nabla²(u) = [0,f] in Omega = ([-1,1] x [-1,1])
%              u(x,y) = g     at d_Omega
%
% Variational problem:
%   a(u,v) = \int_Omega \grad(u) * \sigma \grad(v)
%   f(v)   = \int_Omega f v
%
% Galerkin discretization:
%   \sum_i^num(DOF) u_i \int_Omega \grad(u_i) * \sigma(x,y) \grad(v_j) d(x,y)
% forall
%   j = 1:num(DOF)

%% Set up script.

% Clean up and set verbosity.
clean();
warning('on');
debuging = pick(1, false, true);
verbosity = pick(2, false, true);

%% Set up disctrete Laplace fwd problem.

if debuging
    profile on
end
if verbosity
   fprintf('2D Laplace/Poisson problem\n');
end

% Define (geometry) scaling factor.
scale = pick(1, 1, 50);

% Define outermost grid boundaries.
x = scale * [-1, 1];
y = scale * [-1, 1];

% Define background conductivity.
param_info.val = 1;
param_info.name = {'entire'};

% Define observation points.
RX = pick(1, ...
    [linspace(-0.9, 0.9, scale * 15).', ...
    linspace(0.9, -0.9, scale * 15).'], ... % diagonal profile
    [zeros(scale * 15, 1), ...
    linspace(-0.9, 0.9, scale * 15).']);    % axis parallel profile
RX = scale * RX;

% Define source point and strength.
[TXp, TXd, TXq, TXh] = deal(struct());
TXp.type = 'point_exact';
TXp.coo = scale * pick(1, [0, -1], [0, 0]);
% Note: source AT bnd only reasonable with h. N-BC
TXp.val = 1;                  % discrete    Poisson problem (pole)
TXd.type = 'point_exact';
TXd.coo = scale * pick(1, [.8, .95; -0.02, -0.02], [-0.5, -0.5; 0.5,  0.5]);
TXd.val = [-1, 1];            % discrete    Poisson problem (dipole)
TXq.type = 'point_approx';
TXq.coo = scale * [0.5, 0.5; -0.5, -0.5; -0.5, 0.5; 0.5, -0.5];
TXq.val = [1, 0.5, -1, -0.5]; % discrete    Poisson problem (quadrupole)
TXq.ref_sol = RefSol.getPeak(TXq.val, TXq.coo, 1e-5);
TXh.type = 'reference';
TXh.val = 1;                  % homogeneous Poisson problem
TXh.ref_sol = RefSol.getConst(TXh.val);
%              1    2    3    4
TX = pick(4, TXp, TXd, TXq, TXh);

% Define boundary conditions.
% Note: Detailed preparation follows after setting up the FE system.
[bnd_N, bnd_D, bnd_mix] = deal(struct());
%
bnd_N.type = {'neumann'};
%           xmin xmax ymin ymax
bnd_N.val = {{ 0;   0;   0;   0}};
bnd_N.name = {'xmin', 'xmax', 'ymin', 'ymax'};
bnd_N.quad_ord = 1;
%
bnd_D.type = {'dirichlet'};
%                    xmin xmax ymin ymax
bnd_D.val = {pick(1, {  0;   0;   0;   0}, ... %   homogeneous DRB
                     { 10;  10;   3;   3}, ... % inhomogeneous DRB
                     {  0;   0;  10;   0})};   % inhomogeneous DRB
bnd_D.name = {'xmin', 'xmax', 'ymin', 'ymax'};
%
bnd_mix.type = {'dirichlet', 'neumann'};
%               xmin xmax ymin ymax
bnd_mix.val = {{   0;  0;  [];   0}, ...   % 1 for Dirichlet
               {  [];  [];  0; []}}; ...  % 2 for Neumann
bnd_mix.name = {'xmin', 'xmax', 'ymin', 'ymax'};
bnd_mix.quad_ord = 1;
%                 1      2        3
bnd = pick(3, bnd_N, bnd_D, bnd_mix);

% Choose basic grid type.
mesh_type = 'cube';

% Set number of grid refinements.
ref_steps = 2;

% Set up order of Lagrange elements.
order = pick(2, 1, 2);

% Print status.
if verbosity
   fprintf(sprintf('- use "%s" basic mesh\n', mesh_type));
   fprintf(sprintf('- use "%d" mesh refinements\n', ref_steps));
   fprintf(sprintf('- use "%s" source\n', TX.type));
   if length(bnd.type) > 1
       fprintf('- use mixed boundary conditions\n');
   else
       fprintf(sprintf('- use "%s" boundary conditions\n', bnd.type{1}));
   end
   fprintf(sprintf('- use oder "%d" Lagrange elements\n', order));
end
if verbosity
   fprintf('... FE-FWP set up.\n \n');
end

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, 'bnd', [x, y], ...
    'ref', ref_steps, 'verbosity', verbosity);

%% Set up parameter anomaly.

param = Param.initParam(mesh, param_info);

% Set disturbed area (equals vertical dike).
x_dist = scale * [0, 1];
y_dist = scale * [0.25, 0.55];

% Define parameter of cunductivity (conductor within resistor).
dist = pick(2, param_info.val, 1e3);

% Find all cell midpoints using barycentric coordinates.
lambda_mid = 1/3 + zeros(3, 1);
cell_mid = cellfun(@(x) ...
    {[lambda_mid(1)*x(1,1) + lambda_mid(2)*x(2,1) + lambda_mid(3)*x(3,1); ...
    lambda_mid(1)*x(1,2) + lambda_mid(2)*x(2,2) + lambda_mid(3)*x(3,2)]}, ...
    mesh.cell2cord);

% Find cells, belonging to distrubed area.
cell_dist = cellfun(@(x) (x(1) >= x_dist(1) && x(1) <= x_dist(2)) && ...
        (x(2) >= y_dist(1) && x(2) <= y_dist(2)), cell_mid);

% Set parameter for disturbed area.
param(cell_dist) = dist;

%% Set up FE structure.

fe = FeL.initFiniteElement(order, mesh, RX, verbosity);

%% Set up BC.

bnd = FeL.assignBC(bnd, fe, mesh, param);

%% Set up FEM linear System.

% Set up system matrix.
% (for Poisson/Laplace, this only comprises the stiffness matrix)
sol.A = FeL.assembleStiff(fe, mesh, param, verbosity);

% Set up rhs vector.
sol.b = FeL.assembleRHS(fe, mesh, TX, verbosity);

% Handle boundary conditions.
sol = FeL.treatBC(fe, mesh, sol, bnd, verbosity);
if verbosity
   fprintf('... Linear system and BC set up.\n \n');
end

%% Solve fwd problem.

% Get solution at DOF.
u = FeL.solveFwd(sol, fe, verbosity);

%% Plot solution.

% Solution field.
% Note: in order to obtain a multipol solution, the solution vectors are
% just summed up.
Plot.plotSolution(fe, mesh, sum(u, 2), ...
    'param', param, 'verbosity', verbosity, 'style', '3D');

% Add profile.
% Get solution at RX positions.
phi = fe.I * sum(u, 2);
hold on
    plot3(RX(:,1), RX(:,2), phi, 'r', 'LineWidth', 2);
hold off
if debuging
    profile viewer
end
