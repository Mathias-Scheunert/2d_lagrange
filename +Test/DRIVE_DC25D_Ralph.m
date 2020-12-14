% 2.5D DC problem for a resistive vertical dike in homogenous halfspace.
%
% Problem was introduced by Ralph-Uwe Börner within his lecture to
% visualize the sometimes misleading behavior of the apparent resistivity.
% In this case it indicates a conductive layer with in a resistive
% halfspace whereas a resistive dike in a conductive halfespace is
% considere.
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
%   -\div(\sigma\grad(\phi')) + k_z^2 \sigma \phi' = I \dirac(x_0)
%                    \phi' = 0             at d_Omega_1 (earth interior)
%            d_\phi' / d_n = 0             at d_Omega_2 (surface)
%
% such that
%    \phi = 1/pi \int_{0}^{\inf} \phi' dk_z
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
%
%   To calculate the apparent resistivities, the reciprocity priciple for
%   reducing the number of sources is utilized.
%   -> resp. parts are marked with "% RECI:"

%% Set up script.

clean();

% Set up .geo file
mn1 = 0.02;
mn2 = 0.3;
ab_half = logspace(0.0, 2.0, 17);
Test.createGmsh_Ralph(mn1, mn2, ab_half);

% Set verbosity.
verbosity = pick(2, false, true);

% Define number of uniform grid refinements.
refinement = 0;

% Set order of Lagrange elements.
FE_order = pick(2, 1, 2);

% Define type of numerical integration approach.
FT_type = 'Boerner';

%% Construct Mesh.

% Run Gmsh.
gmsh_path = dir('**/gmsh');
system([gmsh_path.folder, '/gmsh -2 ', ...
        'test.geo -v 0 -format msh2']);
% Load .msh file.
mesh = Mesh.initMesh('gmsh_load', ...
                     'name', 'test.msh', ...
                     'verbosity', verbosity, ...
                     'ref', refinement);
mesh.type = 'gmsh_create';

% Remove obsolete files.
delete('test.geo', 'test.msh');

%% Extract TX/RX.

% Find maximum numbers of TX and RX.
point_idx = mesh.point(mesh.point ~= 0);
str = cellfun(@(x) {x(1)}, mesh.point_names(point_idx));
num = cellfun(@(x) str2double(x(2:end)), mesh.point_names(point_idx));
% RECI: Change roles.
TX_max_num = max(num(strcmp(str, 'M') | strcmp(str, 'N')));
RX_max_num = max(num(strcmp(str, 'A') | strcmp(str, 'B')));

% Construct TX.coo matrix by looping over each pair of A and B and store
% them in an ongoing numbering.
TX = struct();
TX.coo = zeros(TX_max_num*2, 2);
[idx_TX_A, idx_TX_B] = deal(zeros(TX_max_num, 1));
for ii = 1:(TX_max_num)
    % For each TX pair: get current indices of A and B.
% RECI: Change roles.
    cur_A = find(strcmp(mesh.point_names(point_idx), sprintf('M%d', ii)));
    cur_B = find(strcmp(mesh.point_names(point_idx), sprintf('N%d', ii)));
    % Construct index to fill in matrix.
    idx = (ii*2) - 1;
    TX.coo(idx,:) = mesh.vertices(mesh.point2vtx(cur_A),:);
    TX.coo(idx + 1,:) = mesh.vertices(mesh.point2vtx(cur_B),:);
    idx_TX_A(ii) = idx;
    idx_TX_B(ii) = idx + 1;
end
TX.type = 'point_exact';
TX.val = ones(size(TX.coo, 1), 1);

% Construct TX.coo matrix by looping over each pair of M and N and store
% them in an ongoing numbering.
RX = struct();
RX.coo = zeros(RX_max_num*2, 2);
[idx_RX_M, idx_RX_N] = deal(zeros(RX_max_num, 1));
for jj = 1:(RX_max_num)
    % For each RX pair: get current indices of M and N.
% RECI: Change roles.
    cur_M = find(strcmp(mesh.point_names(point_idx), sprintf('A%d', jj)));
    cur_N = find(strcmp(mesh.point_names(point_idx), sprintf('B%d', jj)));
    % Construct index to fill in matrix.
    idx = (jj*2) - 1;
    RX.coo(idx,:) = mesh.vertices(mesh.point2vtx(cur_M),:);
    RX.coo(idx + 1,:) = mesh.vertices(mesh.point2vtx(cur_N),:);
    idx_RX_M(jj) = idx;
    idx_RX_N(jj) = idx + 1;
end

%% Set up boundary conditions.

% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax
bnd.val = {{[];   0;  0;   0; []; []}, ...   % 1 for Dirichlet
           {0;  [];   [];  []; 0;  0}}; ... % 2 for Neumann
bnd.name = {'ymin', 'ymax', 'xmin', 'xmax', 'dike_xmin', 'dike_xmax'};
bnd.quad_ord = 1;

% Define background conductivity.
param_info.val = 1./[100, 1000, 100];
param_info.name = {'dom1', 'dom2', 'dom3'};

%% Set up parameter vector.

param = Param.initParam(mesh, param_info);

%% Assemble 2.5D DC problem.

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;

[fe, sol, FT_info] = App_DC.assembleDC25D(mesh, param, fwd_params, verbosity);

%% Solve 2.5D DC problem.

u_FE = App_DC.solveDC25D(fe, sol, FT_info, verbosity);

%% Construct observed data.

% Get block solution at RX positions.
u = fe.I * u_FE;

% Loop over receivers.
[rhoa_mn1, rhoa_mn2] = deal(zeros(RX_max_num, 1));
for ii = 1:RX_max_num
    % Calculate k for Schlumberger.
    k1 = (pi/mn1 * (ab_half(ii).^2 - mn1^2/4)).';
    k2 = (pi/mn2 * (ab_half(ii).^2 - mn2^2/4)).';

    % Combine solutions (Dipole TX)
    u_MN1 = u(:,idx_TX_A(1)) - u(:,idx_TX_B(1));
    u_MN2 = u(:,idx_TX_A(2)) - u(:,idx_TX_B(2));

    % Combine solutions (Dipole RX)
    u_AB_MN1 = u_MN1(idx_RX_M(ii)) - u_MN1(idx_RX_N(ii));
    u_AB_MN2 = u_MN2(idx_RX_M(ii)) - u_MN2(idx_RX_N(ii));

    % Calculate rhoa vectors
    rhoa_mn1(ii) = k1 .* u_AB_MN1;
    rhoa_mn2(ii) = k2 .* u_AB_MN2;
end

%% Plot mesh.

Plot.plotMesh(mesh, param);
xlim([min(RX.coo(:,1)) + 50, max(RX.coo(:,1)) - 50]);
ylim([0, 15]);
hold on
plot(RX.coo(:,1), RX.coo(:,2), 'vr', 'LineWidth', 3);
plot(TX.coo(:,1), TX.coo(:,2), 'vb', 'LineWidth', 3);
hold off
title(sprintf('Mesh from Gmsh with %d cells.', fe.sizes.cell));

%% Plot results.

% Analytic rhoas from Ralphs Julia script.
rhoa1_ref = [112.20382538046613
             109.24543623927445
             106.05154234236572
             101.59382907825812
              92.96343969355637
              70.68503180427184
              32.79434116131754
              34.18201687250373
              36.36247313486157
              39.340970216329154
              43.072715722155955
              47.458398639527935
              52.34970068705266
              57.56252717951076
              62.895030037572184
              68.1476405540237
              73.14225419352476];
rhoa2_ref = [288.0704270407296
             277.5148208538289
             267.63312610200313
             255.28756902456078
             232.7311161293497
             175.45462257096963
              78.31979165275362
              82.70278721388698
              88.86370126380355
              96.8562016849646
             106.60152381647744
             117.87842009320431
             130.33818362454997
             143.53870887534524
             156.99022763963978
             170.2058651405659
             182.7499742127384];

% Compare rhoas.
figure();
subplot(2, 1, 1)
    loglog(ab_half.', rhoa_mn1, '.-b', ...
           ab_half.', rhoa_mn2, '.-r', ...
           ab_half.', rhoa1_ref, 'o-b', ...
           ab_half.', rhoa2_ref, 'o-r', ...
           'MarkerSize', 10);
    ylabel('rhoa in Ohm*m');
    legend(sprintf('FE MN=%.0d', mn1), sprintf('FE MN=%.0d', mn2), ...
           sprintf('Ralph MN=%.0d', mn1), sprintf('Ralph MN=%.0d', mn2), ...
           'Location', 'North');
subplot(2, 1, 2)
    rel_err_mn1 = (rhoa1_ref - rhoa_mn1) ./ rhoa1_ref * 100;
    rel_err_mn2 = (rhoa2_ref - rhoa_mn2) ./ rhoa2_ref * 100;
    plot(ab_half.', rel_err_mn1, 'b', ...
         ab_half.', rel_err_mn2, 'r');
    xlabel('AB/2 in m');
    ylabel('Rel. error in %');
    ylim([-5, 5]);
