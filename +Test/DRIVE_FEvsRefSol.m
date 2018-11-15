% Script for testing the FE-Code w.r.t. an analytic test function.
%
% Problem 1 'reference':
%   given: 
%       f(x), \nabla²(f(x))
%   solve:
%   -\nabla²(u) + u = -\nabla²(f(x(DOF)))                   in Omega
%              u(x) = -\nabla²(f(x(bndDOF))) + f(x(bndDOF)) at d_Omega_1
%   optional:
%           d_u/d_n = \grad(f(x))) * n                      at d_Omega_2
%
% Problem 1 'poisson':
%   given: 
%       p(x)
%   solve:
%   -\nabla²(u) = \dirac(x_0)      in Omega
%          u(x) = p(x)             at d_Omega_1
%   optional:
%       d_u/d_n = \grad(p(x))) * n at d_Omega_2

%% Set up script.

% Clean up, set verbosity and script.
clean();
warning('on');
debuging = pick(1, false, true);
verbosity = pick(2, false, true);
plotting = pick(2, false, true);

%% Set up disctrete Laplace fwd problem.

if debuging
    profile on
end
if verbosity
   fprintf('Test 2D Lagrange FE using a known analytic solution\n'); 
end

% Set up order of Lagrange elements.
order = pick(2, 1, 2);

% Set number of grid refinements.
ref_steps = 4;

% Define source.
[TXp, TXr] = deal(struct());
%
TXr.type = 'reference';
TXr.ref_sol_u = RefSol.getSin();
TXr.ref_sol.f = @(X, Y) TXr.ref_sol_u.L(X, Y) + TXr.ref_sol_u.f(X, Y);
%
TXp.type = 'point_exact';
TXp.coo = pick(2, [0, 1], [0, 0]);
TXp.val = 1;
TXp.ref_sol_u = RefSol.getPoisson2D(TXp.coo);
%
TX = pick(1, TXr, TXp);

% Define outermost grid boundaries.
switch TX.type
    case 'reference'
        x = [-4, 4];
        y = [-4, 4];
    case 'point_exact'
        x = [-1, 1];
        y = [-1, 1];
end

% Define boundary conditions.
% Note: Detailed preparation follows after setting up the FE system.
bnd_D = struct();
bnd_D.type = {'dirichlet'};
%                       xmin            xmax            ymin            ymax
bnd_D.val = {{TX.ref_sol_u.f; TX.ref_sol_u.f; TX.ref_sol_u.f; TX.ref_sol_u.f}};
bnd_D.name = {'xmin', 'xmax', 'ymin', 'ymax'};
%
bnd_mix = struct();
bnd_mix.type = {'neumann', 'dirichlet'};
%                         xmin            xmax            ymin            ymax                            
bnd_mix.val = {{TX.ref_sol_u.J;             [];             []; TX.ref_sol_u.J}, ...
               {[];             TX.ref_sol_u.f; TX.ref_sol_u.f;             []}};
bnd_mix.name = {'xmin', 'xmax', 'ymin', 'ymax'};
bnd_mix.quad_ord = TX.ref_sol_u.quad_ord;
%
bnd_basic = pick(2, bnd_D, bnd_mix);

% Define observation points.
n_obs = pick(2, 11, 101);
RX = pick(3, ...
    [linspace(x(1), x(end), n_obs).', ...
        linspace(y(end), y(1), n_obs).'], ...               % diagonal profile
    [zeros(n_obs, 1), linspace(y(1), y(end), n_obs).'], ... % axis parallel profile
    [linspace(x(1), x(end), n_obs).', min(y) + zeros(n_obs, 1)], ... % profile on axis
    []);                                                    % none

% Choose basic grid type.
mesh_type = pick(2, 'rhomb', 'cube');

% Print status.
if verbosity
   fprintf(sprintf('- use "%s" basic mesh\n', mesh_type));
   fprintf(sprintf('- use "%d" mesh refinements\n', ref_steps));
   fprintf(sprintf('- use "%s" source\n', TX.type));
   if length(bnd_basic.type) > 1
       fprintf('- use mixed boundary conditions\n');
   else
       fprintf(sprintf('- use "%s" boundary conditions\n', bnd_basic.type{1}));
   end
   fprintf(sprintf('- use oder "%d" Lagrange elements\n', order));
end
if verbosity
   fprintf('... FE-FWP set up.\n \n'); 
end

%% Run fwp.

% Init mesh.
mesh = Mesh.initMesh(mesh_type, 'bnd', [x, y], ...
    'ref', ref_steps, 'verbosity', verbosity);

%% Set up Parameter.

% Set const. background parameter for grid.
param = 1 + zeros(length(mesh.cell2vtx), 1);

%% Set up FE structure.

fe = FeL.initFiniteElement(order, mesh, RX, verbosity);

%% Set up BC.

bnd = bnd_basic;
bnd = FeL.assignBC(bnd, fe, mesh, param);

%% Set up FEM linear System.

% Set up system matrix.
% (for Poisson/Laplace, this only comprises the stiffness matrix)
sol.A = FeL.assembleStiff(fe, mesh, param, verbosity);
if strcmp(TX.type, 'reference')
    sol.A = sol.A + FeL.assembleMass(fe, mesh, param, verbosity);
end

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
    
%% Plot against reference solution.

% Get reference solution at all DOF.
u_ref = arrayfun(@(x, y) TX.ref_sol_u.f(x, y), ...
    fe.DOF_maps.DOF_coo(:, 1), fe.DOF_maps.DOF_coo(:, 2));

if plotting
    Plot.plotSolution(fe, mesh, u_ref, ...
        'param', param, 'verbosity', verbosity, 'style', '3D');
        set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.25, 0.45, 0.5]);
        title('Ref. sol.');
        drawnow();
    % Add profile.
    if ~isempty(RX)
        phi_ref = arrayfun(@(x, y) TX.ref_sol_u.f(x, y), RX(:, 1), RX(:, 2));
        hold on
            plot3(RX(:,1), RX(:,2), phi_ref, 'r', 'LineWidth', 2);
        hold off
        drawnow();
    end
    % Get z-axis limit.
    cur_fig = gcf();
    z_lim = cur_fig.CurrentAxes.ZLim;

    % Get approx. FE-solution.
    Plot.plotSolution(fe, mesh, u, ...
        'param', param, 'verbosity', verbosity, 'style', '3D');
        set(gcf, 'Units', 'normalized', 'Position', [0.5, 0.25, 0.45, 0.5]);
        title('FE sol.');
        zlim(z_lim);
        drawnow();
    % Add profile.
    if ~isempty(RX)
        phi_FE = fe.I * u;
        hold on
            plot3(RX(:,1), RX(:,2), phi_FE, 'r', 'LineWidth', 2);
        hold off
        drawnow();
    end
    
    % Visualize both w.r.t. profile path.
    if ~isempty(RX)
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
            [~, ~, ylim_min] = find(sort(abs(phi_abs)), 1, 'first');
            ylim_min = 10^floor(log10(ylim_min));
            ylim_max = max(abs(phi_abs));
            ylim_max = 10^ceil(log10(ylim_max));
            ylim([ylim_min, ylim_max]);
            ylabel('abs. error');
            xlabel('profile length [m]');
    end
end
if debuging
    profile viewer
end