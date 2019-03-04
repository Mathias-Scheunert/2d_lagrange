% Script for setting up the 2.5D DC problem.
%
% Problem in 3D:
%      x = [x, y, z]
%   \phi = \phi(x)
%           
%   -\div(\sigma\grad(\phi)) = I \dirac(x_0) in Omega
%                       \phi = 0             at d_Omega_1 (earth interior)
%               d_\phi / d_n = 0             at d_Omega_2 (surface)
%
% Problem in 2.5D (at z = 0!):
%       x = [x, y]
%   \phi' = \phi'(x, k_Z)
%
%   -\div(\sigma\grad(\phi')) + k_z^2 \sigma \phi' = I/2 \dirac(x_0)
%                    \phi' = 0             at d_Omega_1 (earth interior)
%            d_\phi' / d_n = 0             at d_Omega_2 (surface)
%
% such that
%    \phi = 2/pi \int_{0}^{\inf} \phi' dk_z
%
% 2.5D Variational problem:
%   a(u,v,k_z) = \int_Omega \grad(\phi') * \sigma \grad(v) + ...
%                k_z^2 \int_Omega \phi' * \sigma v
%   f(v)        = I \int_{Omega} \dirac(x_0) v
%
% Numerical integration over wavenumbers:
%   \phi = \sum_{l = 1}^{N} w_l \phi'(k_{z,l})
%
% Coordinate system (Kartesian):
%  0,0 ------> 
%      |      x
%      |
%      |
%      v
%       y
%
% REMARKS
%   The 2.5D approach uses the bahavior of a known analytic solution
%   (e.g. point source in 3D half-space) to set up the numerical 
%   integration of the seperate 2D solutions (including choice of 
%   wavenumbers and weights) which forms the 2.5D solution.
%   This can only act as an approximation for the solution of an arbitrary 
%   shaped underground!
%   -> I.e. the solutions should be veryfied by comparing them to a full 3D
%   simulation.
%   Furthermore, this implies that only one single point source can be 
%   treated by this approach such that multi-pole arrangements have to be 
%   simulated by adding up solutions for several sources.
%
%   2.5D approach derived from:    Dey A., Morrison H.F.; 1979
%   num. integration derived from: Ralph-Uwe Börner (pers. communication)
%                                  Bing Z.; 1998 (Dissertation)
%                                  Xu S.; 2000

%% Set up script.

% Clean up and set verbosity.
clean();
warning('on');
debugging = pick(1, false, true);
verbosity = pick(2, false, true);

% Define number of uniform grid refinements.
refinement = 2;

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Define type of numerical integration approach.
FT_type = pick(1, 'Boerner', 'Bing', 'Xu');

%% Set up disctrete DC fwd problem.

% Set up domain boundaries.
% (For given TX/RX these may be adapted)
x = [-1000, 1000];
y = [0, 500];
topo_min = -2;
topo_max = 3;

% Define background conductivity.
param_info = struct();
param_info.val = 1/1000;
param_info.name = {'entire'};

% Define source and receiver locations at earth's surface.
TX = struct();
TX.type = 'point_exact';
TX.coo = [0, 0];
TX.val = 1;
%
n_RX = 17;
RX.coo = pick(2, ...
            [linspace(0, 40, n_RX).', ...
            topo_min + (topo_max - topo_min) .* rand(n_RX, 1)], ...
            [linspace(0, 40, n_RX).', ...
            TX.coo(2) + zeros(n_RX, 1)] ...
         );
RX.coo(ismember(RX.coo(:,1), TX.coo(:,1)),:) = [];
RX.coo = round(RX.coo .* 10) ./ 10;

if debugging
    topo = [];
else
    % Add some arbitrary topography.
    topo = pick(2, ...
           [round([linspace(-350, -45, 20).', topo_min + (topo_max - topo_min) .* rand(20,1)], 2); ...
            round([linspace(45, 350, 20).', topo_min + (topo_max - topo_min) .* rand(20,1)], 2)
           ], ...
           []);
end
clear('topo_max', 'topo_min');

% Define mesh.
% Note: TX/RX positions may not be part of the vertex list of 'cube' mesh.
mesh_type = 'gmsh_create';

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd = struct();
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax 
bnd.val = {{[];   0;  0;   0}, ...   % 1 for Dirichlet
           {0;  [];   [];  []}}; ... % 2 for Neumann
bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};
bnd.quad_ord = 1;

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.topo = topo;
fwd_params.bnd = bnd;
fwd_params.dom_bnd = [x, y];
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;
clear('TX', 'RX', 'bnd', 'FT_type', 'FE_order', ...
      'refinement', 'topo', 'x', 'y');

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, 'bnd', fwd_params.dom_bnd, ...
                'ref', fwd_params.ref, 'verbosity', verbosity, ...
                'topo', fwd_params.topo, 'TX', fwd_params.TX.coo, ...
                'RX', fwd_params.RX.coo, 'dom_name', param_info.name{:});

%% Set up parameter vector.

% Set up parameter vector.
param = Param.initParam(mesh, param_info);

% Set disturbed area (equals vertical dike).
x_dist = [15, 22];

% Define parameter of conductivity (conductor within resistor).
if debugging
    sig_anomaly = param_info.val;
else
    sig_anomaly = pick(1, 1/100, param_info.val);

    % Find all cell midpoints using barycentric coordinates.
    lambda_mid = 1/3 + zeros(3, 1);
    cell_mid = cellfun(@(x) ...
        {[lambda_mid(1)*x(1,1) + lambda_mid(2)*x(2,1) + lambda_mid(3)*x(3,1); ...
        lambda_mid(1)*x(1,2) + lambda_mid(2)*x(2,2) + lambda_mid(3)*x(3,2)]}, ...
        mesh.cell2cord);

    % Find cells, belonging to distrubed area.
    cell_dist = cellfun(@(x) ...
        (x(1) > x_dist(1) && x(1) < x_dist(2)), ...
            cell_mid);

    % Set parameter for disturbed area.
    param(cell_dist) = sig_anomaly;
    mesh.parameter_domain(cell_dist) = 2;

    % Update parameter domain vector.
    mesh.parameter_domain_name = [mesh.parameter_domain_name, 'anomaly'];
end

%% Assemble 2.5D DC problem.

[fe, sol, FT_info] = App_DC.assembleDC25D(mesh, param, fwd_params, verbosity);

%% Solve 2.5D DC problem.

u_FE = App_DC.solveDC25D(fe, sol, FT_info, verbosity);

% Obtain potentials at RX positions.
phi_FE = fe.I * u_FE;

%% Compare with analytic point source at top of homogeneous half-space.

% Plot / compare solutions along the profile.
x_plot = sqrt((fwd_params.RX.coo(:,1) - fwd_params.TX.coo(1,1)) .^2 + ...
              (fwd_params.RX.coo(:,2) - fwd_params.TX.coo(1,2)) .^2);

switch fwd_params.FT_type
    case {'Boerner'}
        fig_nun = 2;
    case {'Xu'}
        fig_nun = 20;
    case {'Bing'}
        fig_nun = 200;
end
          
% Get reference solution in 3D.
ref_3D = RefSol.getElectrodeAtHS(1/param_info.val, ...
                                 fwd_params.TX.val, fwd_params.TX.coo);
phi_3D = arrayfun(@(x, y) ref_3D.f(x, y), ...
              fwd_params.RX.coo(:, 1), fwd_params.RX.coo(:, 2));
ref_dike_3D = RefSol.getElectrodeAtVertDike(1/param_info.val, ...
                                       1/sig_anomaly, ...
                                       x_dist(1), x_dist(2)-x_dist(1), ...
                                       fwd_params.TX.coo, ...
                                       fwd_params.TX.val);
phi_dike_3D = arrayfun(@(x, y) ref_dike_3D.f(x, y), ...
              fwd_params.RX.coo(:, 1), fwd_params.RX.coo(:, 2));

% Get asymptotic solution.
% I.e. the numerical integration of the asymptotic Bessel functions
% which map the 2D behavior of the 3D solution for each wavenumbers.
phi_asy = zeros(size(fwd_params.RX.coo, 1), 1);
for i=1:size(fwd_params.RX.coo, 1)
    r = abs(fwd_params.RX.coo(i,1) - fwd_params.TX.coo(1,1));
    switch fwd_params.FT_type
        case {'Boerner', 'Xu'}
%                 RefSol.get_Uk_HR(k, r_TX2RX, I, rho1);
            phi_asy(i) = (2 / pi) * sum(FT_info.w .* ...
                RefSol.get_Uk_HR(FT_info.k, r, 1, 1/param_info.val).');
        case 'Bing'
            % Not known.
            phi_asy = phi_asy * NaN;
    end
end

figure(fig_nun)
subplot(2, 1, 1)
    semilogy(x_plot, phi_3D, 'ok', ...
             x_plot, phi_dike_3D, '+m', ...
             x_plot, phi_asy, 'xb', ...
             x_plot, phi_FE, '-r');
    xlim([x_plot(1), x_plot(end)]);
    ylim([1e0, 1e2]);
    title(sprintf('2.5D half-space solution using %d wavenumbers (%s)', ...
          FT_info.n, fwd_params.FT_type));
    legend('\phi_{HS 3D}', ...
           '\phi_{dike 3D}', ...
           '\phi_{asyHS 3D}', ...
           '\phi_{FE}');
    ylabel('potential');
subplot(2, 1, 2)
    rel_err_FE = (1 - (phi_FE ./ phi_3D)) * 100;
    rel_err_dike_FE = (1 - (phi_FE ./ phi_dike_3D)) * 100;
    rel_err_asy = (1 - (phi_asy ./ phi_3D)) * 100;
    plot(x_plot, rel_err_asy, 'b', ...
         x_plot, rel_err_dike_FE, 'm', ...
         x_plot, rel_err_FE, 'r');
    xlim([x_plot(1), x_plot(end)]);
    ylim([-3, 3]);
    legend('\phi_{HS 3D} vs. \phi_{asyHS 3D}', ...
           '\phi_{dike 3D} vs. \phi_{FE}', ...
           '\phi_{HS 3D} vs. \phi_{FE}');
    ylabel('rel. error');
    xlabel('profile length');

%% Compare some solutions within the wavenumber domain.

if debugging && ~strcmp(FT_info.type, 'Bing')
    % Get asymptotics in 2D.
    phi_ref_2D = cell(FT_info.n, 1);
    for jj = 1:(FT_info.n)
        phi_ref_2D{jj} = arrayfun(@(r) ...
            RefSol.get_Uk_HR(FT_info.k(jj), r, 1, 1/param_info.val).', x_plot);
    end

    % Compare to 2D FE solutions.
    figure(fig_nun+1)
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0, 0, 1, 1])
    all_k_idx = 1:ceil(FT_info.n / 7):FT_info.n;
    for kk = 1:length(all_k_idx)
        subplot(2, length(all_k_idx), kk)
            cur_k_idx = all_k_idx(kk);
            u_cur = FeL.solveFwd(sol{cur_k_idx}, fe);
            phi_FE_2D = fe.I * u_cur;
            semilogy(x_plot, phi_FE_2D, 'r', x_plot, phi_ref_2D{cur_k_idx}, 'ob');
            title(sprintf('k_{%d} = %.1e', cur_k_idx, FT_info.k(cur_k_idx)));
            if kk == 1
                ylabel('$\tilde{\phi}(k_x)$', 'Interpreter', 'latex');
            end
            legend('\phi_{FE}', '\phi_{bessel}');
        subplot(2, length(all_k_idx), length(all_k_idx) + kk)
            rel_err_2D = (1 - (phi_FE_2D ./ phi_ref_2D{cur_k_idx})) * 100;
            plot(x_plot, rel_err_2D);
            if kk == 4
                title('FE-2D vs. analytic 2D solution for variing k_j');
            end
            ylim([-30, 30]);
            if kk == 1
                ylabel('rel. error');
            end
            xlabel('profile length');
    end
end
