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
% clean();
warning('on');
verbosity = pick(2, false, true);

% Define number of uniform grid refinements.
refinement = 1;

% Set boundary constraint.
bc_type = pick(1, 'D_N', 'DtN_N');

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Define type of numerical integration approach.
FT_type = 'Boerner';

file_name = '+Test/zweischicht.msh';

%% Set up disctrete DC fwd problem.

% Set up domain boundaries.
% (For given TX/RX these may be adapted)
x = [-1000, 1000];
y = [0, 500];
topo_min = -2;
topo_max = 3;

% Define background conductivity.
param_info = struct();
param_info.h = 25;
param_info.val = [1/10, 1/1000];
param_info.name = {'schicht_1', 'schicht_2'};

% Define source and receiver locations at earth's surface.
TX = struct();
TX.type = 'point_exact';
TX.coo = [0, 0];
TX.val = 1;
%
n_RX = 7;
RX = struct();
RX.coo = [logspace(0, log10(100), n_RX).', zeros(n_RX, 1)];

% Define mesh.
% Note: TX/RX positions may not be part of the vertex list of 'cube' mesh.
mesh_type = 'gmsh_load';

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
switch bc_type
    case 'D_N'
        bnd = struct();
        bnd.type = {'dirichlet', 'neumann'};
        %         ymin ymax  xmin xmax 
        bnd.val = {{[];   0;  0;   0}, ...   % 1 for Dirichlet
                   {0;  [];   [];  []}}; ... % 2 for Neumann
        bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};
        bnd.quad_ord = 1;
        
    case 'DtN_N'
        bnd.type = {'dtn', 'neumann'};
        %         ymin       ymax       xmin       xmax 
        dtn_fun = App_DC.getD2NBCVals(TX.coo);
        bnd.val = {{[]; dtn_fun.f; dtn_fun.f; dtn_fun.f}, ...  % 1 for D-t-N operator
                   { 0;        [];        [];        []}};     % 2 for Neumann
        bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};
        bnd.quad_ord = dtn_fun.quad_ord;
end

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;

%% Set up mesh.

mesh = Mesh.initMesh(mesh_type, 'name', file_name, ...
                     'ref', fwd_params.ref, 'verbosity', verbosity);

%% Set up parameter vector.

% Set up parameter vector.
param = Param.initParam(mesh, param_info);

%% Assemble 2.5D DC problem.

[fe, sol, FT_info] = App_DC.assembleDC25D(mesh, param, fwd_params, verbosity);

%% Solve 2.5D DC problem.

u_FE = App_DC.solveDC25D(fe, sol, FT_info, verbosity);

% Obtain potentials at RX positions.
phi_FE = fe.I * u_FE;

%% Compare with analytic point source at top of homogeneous half-space.

switch FT_type
    case 'Boerner'
        fig_num = 1;
    case 'Xu'
        fig_num = 10;
    case 'Bing'
        fig_num = 100;
end

% Plot / compare solutions along the profile.
x_plot = sqrt((fwd_params.RX.coo(:,1) - fwd_params.TX.coo(1,1)) .^2 + ...
              (fwd_params.RX.coo(:,2) - fwd_params.TX.coo(1,2)) .^2);

% Get reference solution in 3D.
ref_3D = RefSol.getElectrodeAt2L(1/param_info.val(1), ...
                                 1/param_info.val(2), param_info.h, ...
                                 fwd_params.TX.val, fwd_params.TX.coo);
phi_3D = arrayfun(@(x, y) ref_3D.f(x, y), ...
              fwd_params.RX.coo(:, 1), fwd_params.RX.coo(:, 2));

% Get asymptotic solution.
% I.e. the numerical integration of the asymptotic Bessel functions
% which map the 2D behavior of the 3D solution for each wavenumbers.
phi_asy = zeros(size(fwd_params.RX.coo, 1), 1);
for i=1:size(fwd_params.RX.coo, 1)
    r = abs(fwd_params.RX.coo(i,1) - fwd_params.TX.coo(1,1));
    switch fwd_params.FT_type
        case {'Boerner', 'Xu'}
            phi_asy(i) = (2 / pi) * sum(FT_info.w .* ...
                          RefSol.get_Uk_2L(FT_info.k, r, ...
                            param_info.h, fwd_params.TX.val, ...
                            1/param_info.val(1), 1/param_info.val(2)).');
        case 'Bing'
            % Not known.
            phi_asy = phi_asy * NaN;
    end
end

figure(fig_num)
subplot(2, 1, 1)
    semilogy(x_plot, phi_3D, 'ok', ...
         x_plot, phi_asy, 'xb', ...
         x_plot, phi_FE, '-r');
    xlim([x_plot(1), x_plot(end)]);
    ylim([5e-2, 5e0]);
    title(sprintf('2.5D 2-layer solution using %d wavenumbers (%s)', ...
          FT_info.n, fwd_params.FT_type));
    legend('\phi_{3D}', ...
           '\phi_{asy}', ...
           '\phi_{FE}');
    ylabel('potential');
subplot(2, 1, 2)
    rel_err_FE = (1 - (phi_FE ./ phi_3D)) * 100;
    rel_err_asy = (1 - (phi_asy ./ phi_3D)) * 100;
    plot(x_plot, rel_err_asy, 'b', ...
         x_plot, rel_err_FE, 'r');
    ylim([-10, 10]);
    xlim([x_plot(1), x_plot(end)]);
    legend('\phi_{ref} vs. \phi_{asy}', ...
           '\phi_{ref} vs. \phi_{FE}');
    ylabel('rel. error');
    xlabel('profile length');

%% Compare some solutions within the wavenumber domain.
return
% Get 2-layer asymptotics in 2D.
phi_ref_2D = RefSol.get_Uk_2L(FT_info.k, x_plot, param_info.h, 1, ...
                       1/param_info.val(1), 1/param_info.val(2));

% Compare to 2D FE solutions.
figure(fig_num+1)
all_k_idx = 1:ceil(FT_info.n / 7):FT_info.n;
for kk = 1:length(all_k_idx)
    subplot(2, length(all_k_idx), kk)
        cur_k_idx = all_k_idx(kk);
        u_cur = FeL.solveFwd(sol{cur_k_idx}, fe);
        phi_FE_2D = fe.I * u_cur;
        semilogy(x_plot, phi_FE_2D, 'r', x_plot, phi_ref_2D(:,cur_k_idx), 'ob');
        title(sprintf('k_{%d} = %.3e', cur_k_idx, FT_info.k(cur_k_idx)));
        xlim([x_plot(1), x_plot(end)]);
        if kk == 1
            ylabel('$\tilde{\phi}(k_x)$', 'Interpreter', 'latex');
        end
        legend('\phi_{FE}', '\phi_{bessel}');
    subplot(2, length(all_k_idx), length(all_k_idx) + kk)
        rel_err_2D = (1 - (phi_FE_2D ./ phi_ref_2D(:,cur_k_idx))) * 100;
        plot(x_plot, rel_err_2D);
        if kk == 4
            title('FE-2D vs. analytic 2D solution for variing k_j');
        end
        xlim([x_plot(1), x_plot(end)]);
        ylim([-20, 20]);
        if kk == 1
            ylabel('rel. error');
        end
        xlabel('profile length');
end