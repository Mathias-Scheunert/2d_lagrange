% Script for setting up the 2.5D DC problem.
%
% Coordinate system:
%  0,0 ------> 
%      |      x
%      |
%      |
%      v
%       y

%% Set up script.

% Clean up and set verbosity.
kill();
warning('on');
debugging = pick(2, false, true);

%% Set up disctrete DC fwd problem.

% Set outermost grid boundaries.
x = [-50, 50];
y = [0, 50];

% Define source and receiver locations.
TX.type = 'point_exact';
TX.coo = [-7.5, 0; ...
           7.5, 0];
TX.val = [-1, 1];
%
RX = linspace(-40, 40, 17);

%% Set up FEM.

% Define number of grid refinements.
refinement = pick(2, 2, 4);

% Set order of Lagrange elements.
order = pick(2, 1, 2);

%% Set up mesh.

mesh = Mesh.initMesh('cube', [x, y], refinement, true);

%% Set up Parameter.

% Set disturbed area (equals vertical dike).
x_dist = [-15, -8];
y_dist = [5, 25];

% Define parameter of cunductivity (conductor within resistor).
sig_background = 1/1000;
sig_anomaly = 1/50;

% Find all cell midpoints using barycentric coordinates.
lambda_mid = 1/3 + zeros(3, 1);
cell_mid = cellfun(@(x) ...
    {[lambda_mid(1)*x(1,1) + lambda_mid(2)*x(2,1) + lambda_mid(3)*x(3,1); ...
    lambda_mid(1)*x(1,2) + lambda_mid(2)*x(2,2) + lambda_mid(3)*x(3,2)]}, ...
    mesh.cell2cord);

% Find cells, belonging to distrubed area.
cell_dist = cellfun(@(x) ...
    (x(1) > x_dist(1) && x(1) < x_dist(2)) && ...
    (x(2) > y_dist(1) && x(2) < y_dist(2)), ...
        cell_mid);
    
% Set background parameter for grid.
param = sig_background + zeros(length(mesh.cell2vtx), 1);

% Set parameter for disturbed area.
param(cell_dist) = sig_anomaly;

% Plot grid.
if debugging
    Plot.plotMesh(mesh, param);
    clb = colorbar;
    clb.Label.String = 'parameter value';
end

%% Set up FE structure.


%% Set up FEM linear System.


%% Solve fwd problem.
