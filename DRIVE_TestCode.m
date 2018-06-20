% Script for testing the FE-Code w.r.t. an analytic test function.
%
% Problem:
%   given: 
%       f(x), nabla²(f(x))
%   solve:
%   -nabla²(u) = -nabla²(f(x(DOF)))    in Omega
%       u(x,y) = -nabla²(f(x(bndDOF))) at d_Omega

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

% Clean up, set verbosity and script.
kill();
warning('on');
debug = pick(1, false, true);
verbosity = pick(2, false, true);
plotting = pick(2, false, true);
convergence = pick(1, false, true); % iterate a sequence of refinements
if convergence
    [debug, verbosity, plotting] = deal(false);
end

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
        linspace(y(end), y(1), n_obs).'], ...               % diagonal profile
    [zeros(n_obs, 1), linspace(x(1), x(end), n_obs).'], ... % axis parallel profile
    []);                                                    % none
if convergence
    RX = [];
end

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
if convergence
    mesh_type = 'cube';
end

% Set number of grid refinements.
if convergence
    ref_steps = 0:5;
    [err_L2, err_num_DOF] = deal(zeros(length(ref_steps), 1));
else
    ref_steps = 2;
end

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

% Iterate (if required) over different uniform refinement steps.
for cur_ref = ref_steps
    
    mesh = Mesh.initMesh(mesh_type, [x, y], cur_ref, verbosity);

    %% Set up Parameter.

    % Set const. background parameter for grid.
    param = 1 + zeros(length(mesh.cell2vtx), 1);

    %% Set up FE structure.

    fe = Fe.initFiniteElement(order, mesh, RX, verbosity);

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

    % Get reference solution at all DOF.
    u_ref = arrayfun(@(x, y) ref_sol.f(x, y), ...
        fe.DOF_maps.DOF_coo(:, 1), fe.DOF_maps.DOF_coo(:, 2));

    if convergence
        if cur_ref == 0
            fprintf(sprintf('Num DOF \t L2 error \n'));
        end
        err_L2(cur_ref + 1) = norm(u - u_ref, 2);
        err_num_DOF(cur_ref + 1) = fe.sizes.DOF;  
        fprintf(sprintf('%d \t %d \n', fe.sizes.DOF, norm(u - u_ref, 2)));
    end
end

if convergence
   figure(1);
   loglog(err_num_DOF, err_L2/err_L2(1), 'x-', ...
       err_num_DOF, flipud(err_num_DOF).^fe.order/err_num_DOF(end - 1).^fe.order);
   legend('||u_{FE} - u_{ref}||_2', sprintf('O(h^%d)', fe.order));
end

if plotting
    % Plot reference solution.
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
        [~, ~, ylim_min] = find(sort(abs(phi_abs)), 1, 'first');
        ylim_min = 10^floor(log10(ylim_min));
        ylim_max = max(abs(phi_abs));
        ylim_max = 10^ceil(log10(ylim_max));
        ylim([ylim_min, ylim_max]);
        ylabel('abs. error');
        xlabel('profile length [m]');
end

if debug
    profile viewer
end