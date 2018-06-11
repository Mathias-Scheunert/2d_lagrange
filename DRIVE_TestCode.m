% Script for testing the FE-Code w.r.t. an analytic test function.
%
% Problem:
%   given: 
%       f(x), nabla²(f(x)
%   solve:
%   -nabla²(u) = -nabla²(f(x(DOF)))    in Omega
%       u(x,y) = -nabla²(f(x(bndDOF))) at d_Omega

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
   fprintf('Test 2D Lagrange FE using a known analytic solution\n'); 
end

% Define outermost grid boundaries.
x = [-4, 4];
y = [-4, 4];

% Define observation points.
n_obs = pick(2, 11, 101);
RX = pick(1, ...
    [linspace(x(1), x(end), n_obs).', ...
        linspace(y(end), y(1), n_obs).'], ...            % diagonal profile
    [zeros(n_obs, 1), linspace(x(1), x(end), n_obs).']); % axis parallel profile

% Define source.
TX.type = 'reference';
ref_sol = RefSol.getSinFunction();
TX.ref_sol.f = ref_sol.L;

% Define boundary conditions.
bnd = struct();
bnd.type = 'dirichlet';
bnd.val = [];

% Choose basic grid type.
mesh_type = pick(2, 'rhomb', 'cube');

% Set number of grid refinements.
ref_steps = 6;

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
   fprintf(sprintf('- use oder "%d"-order Lagrange elements\n', order));
end
if verbosity
   fprintf('... FE-FWP set up.\n \n'); 
end

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, [x, y], ref_steps, verbosity);
if verbosity
   fprintf('... Mesh struct initialized.\n \n'); 
end

%% Set up Parameter.
 
% Set const. background parameter for grid.
param = 1 + zeros(length(mesh.cell2vtx), 1);

%% Set up FE structure.

fe = Fe.initFiniteElement(order, mesh, RX, verbosity);
if verbosity
   fprintf('... FE struct initialized.\n \n'); 
end

%% Set up reference solution quantities.

% Treat inhomogeneous Dirichlet boundary values.
% Get DOF and it's coords at bnd.
bnd.bndDOF = Fe.getBndDOF(fe, mesh);
% Get desired reference function values at bnd.
bnd.val = arrayfun(@(x, y) ref_sol.f(x, y), ...
    bnd.bndDOF.bnd_DOF_coo(:, 1), bnd.bndDOF.bnd_DOF_coo(:, 2));

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
[sol, bnd] = Fe.treatDirichlet(fe, mesh, sol, bnd);
if verbosity
   fprintf('done.\n'); 
end
if verbosity
   fprintf('... Linear system and BC set up.\n \n'); 
end

%% Solve fwd problem.

% Get solution at DOF.
u = Fe.solveFwd(sol, fe, verbosity);
if verbosity
   fprintf('... FWP solved.\n \n'); 
end

%% Plot solution.

% Plot reference solution.
u_ref = arrayfun(@(x, y) ref_sol.f(x, y), ...
    fe.DOF_maps.DOF_coo(:, 1), fe.DOF_maps.DOF_coo(:, 2));
Plot.plotSolution(fe, mesh, u_ref, param, verbosity);
    set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.25, 0.45, 0.5]);
    title('Ref. sol.');
    drawnow();
% Add profile.
phi_ref = arrayfun(@(x, y) ref_sol.f(x, y), RX(:, 1), RX(:, 2));
hold on
    plot3(RX(:,1), RX(:,2), phi_ref, 'r', 'LineWidth', 2);
hold off
drawnow();

% Get approx. FE-solution.
Plot.plotSolution(fe, mesh, u, param, verbosity);
    set(gcf, 'Units', 'normalized', 'Position', [0.5, 0.25, 0.45, 0.5]);
    title('FE sol.');
    drawnow();
% Add profile.
phi_FE = fe.I * u;
hold on
    plot3(RX(:,1), RX(:,2), phi_FE, 'r', 'LineWidth', 2);
hold off
drawnow();


% Visualize both w.r.t. profile path.
x = sqrt((RX(:,1) - RX(1,1)) .^2 + (RX(:,2) - RX(1,2)) .^2);
figure();
set(gcf, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
subplot(2, 1, 1)
    plot(x, phi_FE, '-or', ...
        x, phi_ref, '--xk');
    xlim([0, x(end)]);
    ylabel('sol. value');
    legend('FE-sol', 'Ref-sol', 'Location', 'best');
% Show relative errors.
phi_abs = (phi_FE - phi_ref);
subplot(2, 1, 2)
    semilogy(x, abs(phi_abs));
    xlim([0, x(end)]);
    ylim([1e-8, 1e0]);
    ylabel('abs. error');
    xlabel('profile length [m]');
    
if verbosity
   fprintf('... Solution visualized.\n'); 
end
if debug
    profile viewer
end