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
verbosity = pick(2, false, true);

% Enable small parameter anomaly.
include_anomaly = pick(1, false, true);

% Enable smooth topography.
include_topo = pick(1, false, true);

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

if include_topo
    % Define some smooth topography model.
    elev_x = [-100, -90, -80, -65,  -30, -10,   0, 15, 35, 65, 85, 100];
    elev_y = [   3,   4,   3,  -1, -2.5,  -8, -10, -9, -7, -5, -8, -10];
    warning('off');
    topo_poly = polyfit(elev_x, elev_y, length(elev_x) - 3);
    warning('on');

    % Derive topography info.
    topo_x = -100:100;
    topo_y = topo_poly(end);
    for ii = 1:length(topo_poly) - 1
        topo_y = topo_y + (topo_poly(ii) * topo_x.^(length(topo_poly) - ii));
    end
    topo = [topo_x.', topo_y.'];
    
else
    % Set electrode positions at top of halfspace.
    topo = [(-100:100).', 0*(-100:100).'];
end

% Define electrode positions.
% (Considering the TUBAF multi electrode system with 121 electrodes,
% assuming an topography spacing of 1m!)
ixd_ele_in_topo = topo(:,1) > -61 & topo(:,1) < 61;
ele_pos = topo(ixd_ele_in_topo, :);

% Create a measurement configuration.
conf_type = 'wenner';
dc_conf = App_DC.createConfigBERT(ele_pos, conf_type, 'verbosity', verbosity);
clear('elev_x', 'elev_y', 'ele_pos', 'ixd_ele_in_topo', ...
      'topo_poly', 'topo_x', 'topo_y', 'ii');

% Define mesh.
% Note: TX/RX positions may not be part of the vertex list of 'cube' & 
% 'rhomb' mesh.
mesh_type = 'gmsh_create';

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax 
bnd.val = {{[];   0;  0;   0}, ...   % 1 for Dirichlet
           {0;  [];   [];  []}}; ... % 2 for Neumann
bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};
bnd.quad_ord = 1;

% Define background conductivity.
param_info.val = 1/1000;
param_info.name = {'entire'};

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = dc_conf.TX;
fwd_params.RX = dc_conf.RX;
fwd_params.topo = topo;
fwd_params.bnd = bnd;
fwd_params.dom_bnd = [x, y];
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;
clear('TX', 'RX', 'bnd', 'FT_type', 'FE_order', ...
      'refinement', 'topo', 'x', 'y');

%% Set up mesh.

warning('off', 'Mesh:createGmsh:bnd2Close');
mesh = Mesh.initMesh(mesh_type, 'bnd', fwd_params.dom_bnd, ...
    'ref', fwd_params.ref, 'verbosity', verbosity, ...
    'topo', fwd_params.topo, 'TX', fwd_params.TX.coo, ...
    'RX', fwd_params.RX.coo, 'dom_name', param_info.name{:});

%% Set up parameter vector.

% Set up parameter vector.
param = Param.initParam(mesh, param_info);

% Set disturbed area (equals vertical dike).
x_dist = [-80, -10];
y_dist = [10, 80];

% Define parameter of conductivity (conductor within resistor).
if include_anomaly
    sig_anomaly = 10;

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

% Obtain vector of observed data for the 
rhoa = App_DC.getObsData(u_FE, fe, dc_conf);

%% Show some results.

% Visualize Mesh and TX/RX positions.
Plot.plotMesh(mesh, param);
xlim([-110, 110]);
ylim([-15, 30]);
hold on
plot(dc_conf.TX.coo(:,1), dc_conf.TX.coo(:,2), 'vy');
hold off

% Error between apparent and true resistivity for halfspace.
rho_true = 1/param_info.val + 0*rhoa;
rel_err_rho = ((rho_true - rhoa) ./ rho_true) * 100;
figure();
plot(1:length(rel_err_rho), rel_err_rho);
title('Relative error between observed rhoa and true halfspace resistivity');
xlabel('Number of observed data in survey');
ylabel('Error in %');
xlim([0, length(rhoa)]);
ylim([-2, 2]);