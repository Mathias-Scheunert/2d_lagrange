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
refinement = 0;

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Define type of numerical integration approach.
FT_type = 'Boerner';

%% Set up disctrete DC fwd problem.

% Set up domain boundaries.
domain_bnd = [-1000, 1000, 0 ,500];
dike_info = [15, 7];

% Define source and receiver locations at earth's surface.
TX = struct();
TX.type = 'point_exact';
TX.coo = [0, 0];
TX.val = 1;
%
n_RX = 17;
RX = struct();
% RX.coo = [logspace(0, log10(30), n_RX).', zeros(n_RX, 1)];
RX.coo = [linspace(0, 40, n_RX).', zeros(n_RX, 1)];
RX.coo(ismember(RX.coo(:,1), TX.coo(:,1)),:) = [];
RX.coo = round(RX.coo .* 10) ./ 10;

% Define mesh.
% Note: TX/RX positions may not be part of the vertex list of 'cube' mesh.
file_name = 'vert_dike';
Mesh.createVerticalDikeMesh(domain_bnd, [TX.coo; RX.coo], ...
                       dike_info(1), dike_info(2), ...
                       file_name);

% Define background conductivity.
param_info = struct();
param_info.val = [1/1000, 1/100, 1/1000];
param_info.name = {'left', 'middle', 'right'};

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd = struct();
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax 
bnd.val = {{[];   0;  0;   0}, ...   % 1 for Dirichlet
           {0;  [];   [];  []}}; ... % 2 for Neumann
bnd.name = {'top', 'bot', 'left', 'right'};
bnd.quad_ord = 1;

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;

%% Set up mesh.

mesh_type = 'gmsh_load';
mesh = Mesh.initMesh(mesh_type, 'name', [file_name, '.msh'], ...
                                'ref', fwd_params.ref, ...
                                'verbosity', verbosity);

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
        fig_num = 3;
    case 'Xu'
        fig_num = 30;
    case 'Bing'
        fig_num = 300;
end

% Plot / compare solutions along the profile.
x_plot = sqrt((fwd_params.RX.coo(:,1) - fwd_params.TX.coo(1,1)) .^2 + ...
              (fwd_params.RX.coo(:,2) - fwd_params.TX.coo(1,2)) .^2);

% Get reference solutions in 3D.
ref_dike_3D = RefSol.getElectrodeAtVertDike(1/param_info.val(1), ...
                                       1/param_info.val(2), ...
                                       dike_info(1), dike_info(2), ...
                                       fwd_params.TX.coo, ...
                                       fwd_params.TX.val);
phi_dike_3D = arrayfun(@(x, y) ref_dike_3D.f(x, y), ...
              fwd_params.RX.coo(:, 1), fwd_params.RX.coo(:, 2));
ref_HS_3D = RefSol.getElectrodeAtHS(1/param_info.val(1), ...
                                    fwd_params.TX.val, fwd_params.TX.coo);
phi_HS_3D = arrayfun(@(x, y) ref_HS_3D.f(x, y), ...
                  fwd_params.RX.coo(:, 1), fwd_params.RX.coo(:, 2));

figure(fig_num)
subplot(2, 1, 1)
    semilogy(x_plot, phi_dike_3D, 'ok', ...
             x_plot, phi_HS_3D, 'xb', ...
             x_plot, phi_FE, '-r');
    xlim([x_plot(1), x_plot(end)]);
    ylim([1e0, 1e2]);
    title(sprintf('2.5D vertical dike solution using %d wavenumbers (%s)', ...
          FT_info.n, fwd_params.FT_type));
    legend('\phi_{dike 3D}', ...
           '\phi_{HS 3D}', ...
           '\phi_{FE}');
    ylabel('potential');
subplot(2, 1, 2)
    rel_err_FE = (1 - (phi_FE ./ phi_dike_3D)) * 100;
    rel_err_HS = (1 - (phi_HS_3D ./ phi_dike_3D)) * 100;
    plot(x_plot, rel_err_HS, 'b', ...
         x_plot, rel_err_FE, 'r');
    ylim([-10, 10]);
    xlim([x_plot(1), x_plot(end)]);
    legend('\phi_{dike ref} vs. \phi_{HS ref}', ...
           '\phi_{dike ref} vs. \phi_{FE}');
    ylabel('rel. error');
    xlabel('profile length');

return;
%%
Plot.plotMesh(mesh), ylim([-10, 20]), xlim([-20, 60]);
hold on; 
plot(RX.coo(:,1), RX.coo(:,2), '*b');