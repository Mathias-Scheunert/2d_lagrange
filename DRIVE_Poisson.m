% Script for setting up the 2D Poisson problem.
%
% Problem:
%   -\grad(\sigma\grad(u)) = ...
%   -\sigma\nabla²(u) = [0,f] in Omega = ([-1,1] x [-1,1])
%              u(x,y) = g     at d_Omega
%
% Variational problem:
%   a(u,v) = \int_Omega \grad(u) * \sigma \grad(v) d(x,y)
%   f(v)   = \int_Omega f v               d(x,y)
%
% Galerkin discretization:
%   \sum_i^num(DOF) u_i \int_Omega \grad(u_i) * \sigma(x,y) \grad(v_j) d(x,y) 
% forall
%   j = 1:num(DOF)

%% Set up script.

% Clean up and set verbosity.
kill();
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
TXp.coo = scale * pick(2, [0, 1], [0, 0]); 
% Note: source AT bnd only reasonable with h. N-BC
TXp.val = 1;                  % discrete    Poisson problem (pole)
TXd.type = 'point_exact';
TXd.coo = scale * pick(1, [.8, .95; -0.02, -0.02], [-0.5, -0.5; 0.5,  0.5]);             
TXd.val = [-1, 1];            % discrete    Poisson problem (dipole)
TXq.type = 'point_approx';
TXq.coo = scale * [0.5, 0.5; -0.5, -0.5; -0.5, 0.5; 0.5, -0.5];
TXq.val = [1, 0.5, -1, -0.5]; % discrete    Poisson problem (quadrupole)
TXq.ref_sol = RefSol.getPeakFunction(TXq.val, TXq.coo, 1e-5);
TXh.type = 'reference';
TXh.val = 1;                  % homogeneous Poisson problem
TXh.ref_sol = RefSol.getConstFunction(TXh.val);
%              1    2    3    4
TX = pick(1, TXp, TXd, TXq, TXh);

% Choose basic grid type.
mesh_type = pick(2, 'rhomb', 'cube', 'external');

% Define boundary conditions.
% Note: Detailed preparation follows after setting up the FE system.
[bnd_N, bnd_D, bnd_mix] = deal(struct());
%
bnd_N.type = {'neumann'};
% Note:
% x (left -> right)
% y (bottom -> top)
%                   bot top left right
bnd_N.val = {pick(2, {0;  0;  0;  0}, {0;  0;  1;  -1})};
%
bnd_D.type = {'dirichlet'};
%                     bot top left right
bnd_D.val = {pick(1, {  0;  0;   0;    0 }, ... %   homogeneous DRB
                     {  3;  3;  10;   10 }, ... % inhomogeneous DRB
                     { 10;  0;   0;    0 })};   % inhomogeneous DRB
%
bnd_mix.type = {'dirichlet', 'neumann'};
%                       bot top left right
bnd_mix.val = pick(2,{{ 10; [];  3;    [] }, ...  % 1 for Dirichlet
                      { [];  0; [];     0 }}, ... % 1 for Neumann
                     {{ 10; 10;  [];    [] }, ...  % 2 for Dirichlet
                      { [];  []; 1e-2; 1e-2 }}); ...% 2 for Neumann
%                 1      2        3      
bnd = pick(1, bnd_N, bnd_D, bnd_mix);

% Set number of grid refinements.
ref_steps = 0;

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

mesh = Mesh.initMesh(mesh_type, [x, y], ref_steps, verbosity);

%% Set up Parameter.

% TODO: if external mesh is provided, also these info has to be given.

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

%% Set up BC.

bnd = Fe.assignBC(bnd, fe, mesh, param);

%% Set up FEM linear System.

% Set up system matrix.
% (for Poisson/Laplace, this only comprises the stiffness matrix)
sol.A = Fe.assembleStiff(fe, mesh, param, verbosity);

% Set up rhs vector.
sol.b = Fe.assembleRHS(fe, mesh, TX, verbosity);

% Handle boundary conditions.
[sol, bnd] = Fe.treatBC(fe, mesh, sol, bnd, verbosity);
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
if debuging
    profile viewer
end