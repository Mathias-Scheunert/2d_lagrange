% Script to test meshing for 2D.

%% Prepare Script.

% Clean up and set verbosity.
kill();
warning('on');
debug = pick(2, false, true);

%% Define domain boundaries.

% Set outermost grid boundaries.
x = [0, 10];
y = [2, 5];

% Choose basic grid type.
grid_type = pick(1, 'rhomb', 'cube');

% Define number of grid refinements.
refinement = pick(1, 0, 1, 2);

%% Set up grid.

% Construct triangular grid for a rectangular domain.
mesh = Mesh.initMesh(grid_type, [x, y], 0);
n_cells = length(mesh.cell2vtx);

% Plot grid.
if debug
    Plot.plotMesh(mesh, (1:n_cells)/n_cells, false);
    clb = colorbar;
    clb.Label.String = 'cell number';
    clb_labels = cellfun(@(x) {num2str(x)}, num2cell(round(clb.Ticks * n_cells)));
    clb.Ticks = round(clb.Ticks * n_cells) / n_cells;
    clb.TickLabels = clb_labels;
    drawnow();
end

%% Refine grid.

% Apply refinement.
mesh = Mesh.refineMeshUniform(mesh, 1);
n_cells = length(mesh.cell2vtx);

% Plot grid.
if debug
    Plot.plotMesh(mesh, (1:n_cells)/n_cells, true);
    clb = colorbar;
    clb.Label.String = 'cell number';
    clb_labels = cellfun(@(x) {num2str(x)}, num2cell(round(clb.Ticks * n_cells)));
    clb.Ticks = round(clb.Ticks * n_cells) / n_cells;
    clb.TickLabels = clb_labels;
    drawnow();
end

%% Use Affine map to operate on the grid (triangles).

% Set test_point.
%                    1      2      3           4
test_point = pick(4, [2,3], [5,3], [2.5, 3.5], [9.437, 5]);

% Get barycentric coordinates for each triangle w.r.t. test_point.
lambda_all = cellfun(@(x,y) {Mesh.getAffineMap(x, mesh, test_point)}, ...
    num2cell(1:n_cells).');
% Expand to the third barycentric coordinate.
% lambda_all = cellfun(@(x) {[x.xy_ref.'; 1 - x.xy_ref(1) - x.xy_ref(2)]}, lambda_all);
lambda_all = cellfun(@(x) {x.xy_ref.'}, lambda_all);

% Use barycentric coordinates to find cell_num(s) containing test_point:
% Check affiliation: Only if all three lambdas are >= 0, the point belongs 
% to the triangle.
% (note: as only two coordinates are used, the test differs slightly)
cells_fit = cellfun(@(x) all(x >= -eps & x <= 1 + eps) && (sum(x) - 1 < eps), ...
    lambda_all);

% Show point and obtained cells.
if debug
    Plot.plotMesh(mesh, double(cells_fit), false);
    hold on
    plot(test_point(1), test_point(2), '+r', ...
        'MarkerSize', 8, 'LineWidth', 2);
    hold off
end
drawnow();

% Find all cell midpoints using barycentric coordinates.
lambda_mid = 1/3 + zeros(3, 1);
cell_mid = cellfun(@(x) ...
    {[lambda_mid(1)*x(1,1) + lambda_mid(2)*x(2,1) + lambda_mid(3)*x(3,1); ...
    lambda_mid(1)*x(1,2) + lambda_mid(2)*x(2,2) + lambda_mid(3)*x(3,2)]}, ...
    mesh.cell2cord);

% Add mid points to plot.
if debug
    hold on
    cellfun(@(x) plot(x(1), x(2), '*k'), cell_mid);
    hold off
end
drawnow();