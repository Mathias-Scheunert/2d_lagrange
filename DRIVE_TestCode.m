% Script for testing the FE-Code w.r.t. an analytic test function.
%
% Problem:
%   given: 
%       f(x), \nabla²(f(x))
%   solve:
%   -\nabla²(u) = -\nabla²(f(x(DOF)))    in Omega
%        u(x,y) = -\nabla²(f(x(bndDOF))) at d_Omega
% Variants:
%   f(x) = \dirac(x_0)
%   f(x) = -\nabla²(f(x(DOF))) 

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
convergence = pick(2, false, true); % iterate a sequence of refinements
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

% Set up order of Lagrange elements.
if convergence
    order = [1, 2];
else
    order = pick(1, 1, 2);
end

% Define source.
[TXp, TXr] = deal(struct());
%
TXr.type = 'reference';
TXr.ref_sol_u = RefSol.getSinFunction();
TXr.ref_sol.f = TXr.ref_sol_u.L;
%
TXp.type = 'point_exact';
TXp.coo = pick(2, [0, 1], [0, 0]);
TXp.val = 1;
TXp.ref_sol_u.f = @(x, y) -1/(2*pi) * log(norm([x; y] - TXp.coo(:)));
x_sym = sym('x', 'real');
y_sym = sym('y', 'real');
TXp.ref_sol_u.grad = [diff(TXp.ref_sol_u.f, x_sym); diff(TXp.ref_sol_u.f, y_sym)];
TXp.ref_sol_u.grad = matlabFunction(TXp.ref_sol_u.grad, 'Vars', {'x', 'y'});
clear('x_sym', 'y_sym');
TXp.ref_sol_u.quad_ord = 4;
%
TX = pick(2, TXr, TXp);

% Define outermost grid boundaries.
switch TX.type
    case 'reference'
        x = [-4, 4];
        y = [-4, 4];
    case 'point_exact'
        x = [-1, 1];
        y = [-1, 1];
end

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

% Define (inhomogeneous Dirichlet) boundary conditions.
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
    ref_steps = 1:4;
    [err_L2, err_H1, err_num_DOF] = deal(cell(length(order), 1));
else
    ref_steps = 3;
end

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

% Iterate (if required) over different Lagrange orders.
for cur_order = order
    % Preallocate variables.
    [err_L2{cur_order}, err_H1{cur_order}, err_num_DOF{cur_order}] = ...
        deal(zeros(length(ref_steps), 1));
    % Iterate (if required) over different uniform refinement steps.
    for cur_ref = ref_steps
        % Print status.
        if convergence
            fprintf('Test "%d" order with "%d" refinements ...', ...
                cur_order, cur_ref);
        end
        
        % Init mesh.
        mesh = Mesh.initMesh(mesh_type, [x, y], cur_ref, verbosity);

        %% Set up Parameter.

        % Set const. background parameter for grid.
        param = 1 + zeros(length(mesh.cell2vtx), 1);

        %% Set up FE structure.

        fe = Fe.initFiniteElement(cur_order, mesh, RX, verbosity);

        %% Set up reference solution quantities.

        % Treat inhomogeneous Dirichlet boundary values.
        % Get DOF and it's coords at bnd.
        bnd.bndDOF = Fe.getBndDOF(fe, mesh);
        % Get desired reference function values at bnd.
        bnd.val = arrayfun(@(x, y) TX.ref_sol_u.f(x, y), ...
            bnd.bndDOF.bnd_DOF_coo(:, 1), bnd.bndDOF.bnd_DOF_coo(:, 2));

        %% Set up FEM linear System.

        % Set up system matrix.
        % (for Poisson/Laplace, this only comprises the stiffness matrix)
        sol.A = Fe.assembleStiff(fe, param, verbosity);

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

        %% Handle errors.

        if convergence
            cur_idx = cur_ref - ref_steps(1) + 1;
            [err_L2{cur_order}(cur_idx), err_H1{cur_order}(cur_idx)] = ...
                Fe.getError(mesh, fe, u, TX.ref_sol_u);
            err_num_DOF{cur_order}(cur_idx) = fe.sizes.DOF;  
        end
        
        if convergence
            fprintf(' done. \n');
        end
    end
end

%% Plot solution.

if convergence
    for kk = 1:length(order)
        fprintf(sprintf('Order %d: \n Num DOF  L2 error \t H1 error \n', ...
            order(kk)));
        for ii = 1:length(ref_steps)
            fprintf(sprintf('%d \t %d \t %d \n', ...
                err_num_DOF{kk}(ii), err_L2{kk}(ii), err_H1{kk}(ii)));
        end
    end
end
    
% Get reference solution at all DOF.
u_ref = arrayfun(@(x, y) TX.ref_sol_u.f(x, y), ...
    fe.DOF_maps.DOF_coo(:, 1), fe.DOF_maps.DOF_coo(:, 2));

if convergence
   figure(1);
   for kk = 1:length(order)
       subplot(2, 1, kk)
           h_x = logspace(log10(err_num_DOF{kk}(1)), 5, 100);
           h = fliplr(h_x);
           loglog(err_num_DOF{kk}, err_L2{kk}/err_L2{kk}(1), 'x-', ...
               err_num_DOF{kk}, err_H1{kk}/err_H1{kk}(1), 'o-', ...
               h_x, ones(size(h)), '-b', ...
               h_x, h/h(1), '-k' ...
               );
           ylim([1e-3, 1e1]);
           title(sprintf('u_{FE} vs. u_{ref} w.r.t. DOF (%d. order Lagrange)', ...
               order(kk)));
           ylabel('error');
           xlabel('DOF');
           legend('||u_{FE} - u_{ref}||_{L2}', '||u_{FE} - u_{ref}||_{H1}', ...
               'O(const)', 'O(h)', ...
               'Location', 'EastOutside');
   end
end

if plotting
    % Plot reference solution.
    Plot.plotSolution(fe, mesh, u_ref, param, verbosity);
        set(gcf, 'Units', 'normalized', 'Position', [0.05, 0.25, 0.45, 0.5]);
        title('Ref. sol.');
        drawnow();
    % Add profile.
    phi_ref = arrayfun(@(x, y) TX.ref_sol_u.f(x, y), RX(:, 1), RX(:, 2));
    hold on
        plot3(RX(:,1), RX(:,2), phi_ref, 'r', 'LineWidth', 2);
    hold off
    drawnow();
    % Get z-axis limit.
    cur_fig = gcf();
    z_lim = cur_fig.CurrentAxes.ZLim;

    % Get approx. FE-solution.
    Plot.plotSolution(fe, mesh, u, param, verbosity);
        set(gcf, 'Units', 'normalized', 'Position', [0.5, 0.25, 0.45, 0.5]);
        title('FE sol.');
        zlim(z_lim);
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