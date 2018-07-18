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
%   num. integration derived from: Bing Zhou; 1998 (Dissertation)
%                                  Ralph-Uwe Börner (pers. communication)

%% Set up script.

% Clean up and set verbosity.
kill();
warning('on');
debugging = pick(2, false, true);
verbosity = pick(2, false, true);

%% Set up disctrete DC fwd problem.

% Set up domain boundaries.
x = [-50, 50];
y = [0, 50];

% Define source and receiver locations at earth's surface.
TX.type = 'point_exact';
TX.coo = [0, 0];
TX.val = 1;
%
RX.coo = [linspace(-40, 40, 17).', zeros(17, 1)];
RX.coo(ismember(RX.coo, TX.coo, 'rows'),:) = [];

% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
%           xmin xmax ymin ymax
bnd.val = {{   0;   0;  [];   0}, ...   % 1 for Dirichlet
           {  [];  [];   0;  []}}; ...  % 2 for Neumann
bnd.quad_ord = 1;

%% Set up FEM.

% Define number of grid refinements.
refinement = pick(2, 2, 4, 5);

% Set order of Lagrange elements.
order = pick(2, 1, 2);

%% Set up mesh.

mesh = Mesh.initMesh('cube', [x, y], refinement, verbosity);

%% Set up conductivity distribution.

% Set disturbed area (equals vertical dike).
x_dist = [-15, -8];
y_dist = [5, 25];

% Define parameter of conductivity (conductor within resistor).
if debugging
    [sig_anomaly, sig_background] = deal(2 * pi);
else
    sig_background = 1/1000;
    sig_anomaly = 1/50;
end

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
    
% Set background parameter for grid.
param = sig_background + zeros(length(mesh.cell2vtx), 1);

% Set parameter for disturbed area.
param(cell_dist) = sig_anomaly;

%% Set up FE structure.

fe = Fe.initFiniteElement(order, mesh, RX.coo, verbosity);
bnd = Fe.assignBC(bnd, fe, mesh, param);

%% Treat 2.5D wavenumber domain handling and set up DC-FE system.

% Set up invariant rhs vector.
% Note: source strengh needs to be halved.
TX.val = TX.val / 2;
rhs = Fe.assembleRHS(fe, mesh, TX, verbosity);

% Set up invariant system matrix parts.
A_GradDiv = Fe.assembleStiff(fe, mesh, param, verbosity);
A_Mass = Fe.assembleMass(fe, mesh, param, verbosity);

% Define number of wavenumbers for both parts of the asymptotic behavior of
% the analytic solution.
n_k_log = 20;
n_k_exp = 10;
n_wavenum = n_k_log + n_k_exp;

% Get quadrature rules for both asymptotics.
[x_Leg, w_Leg] = Quad.getQuadratureRule(n_k_log, 1, ...
    'type', 'Legendre', 'bnd', [0, 1]);
[x_Lag, w_Lag] = Quad.getQuadratureRule(n_k_exp, 1, ...
    'type', 'Laguerre');

% Get critical wavenumber (transition of asymptotic behavior).
r_TX2RX = sqrt((TX.coo(1) - RX.coo(:,1)).^2 + (TX.coo(2) - RX.coo(:,2)).^2);
k_crit = 1 / (2 * min(r_TX2RX));

% Derive wavenumbers and weights.
k_log = k_crit * x_Leg(:).^2;
w_log = 2 * k_crit * (x_Leg(:) .* w_Leg(:));
k_exp = k_crit * (x_Lag(:) + 1);
w_exp = k_crit * (exp(x_Lag(:)) .* w_Lag(:));
k = [k_log; k_exp];
w = [w_log; w_exp];

% Set up FE systems for all wave numbers.
if verbosity
    fprintf('Assemble linear system and incorporate BC ... '); 
end
sol = cell(n_wavenum, 1);
for ii = 1:(n_wavenum)
    % Set up system matrix
    sol{ii}.A = A_GradDiv + k(ii)^2 * A_Mass;

    % Set up rhs.
    sol{ii}.b = rhs;

    % Handle boundary conditions.
    sol{ii} = Fe.treatBC(fe, mesh, sol{ii}, bnd);
end
if verbosity
    fprintf('done.\n');
    fprintf('... DC 2.5D problem initialized.\n \n');
end
%% Solve fwd problems.

if verbosity
    fprintf('Solve linear systems for all wavenumbers ... '); 
end
u = cell(n_wavenum, 1);
for jj = 1:(n_wavenum)
    u{jj} = Fe.solveFwd(sol{jj}, fe);
end
if verbosity
    fprintf('done.\n'); 
end

%% Apply numerical integration over wavenumber domain.

if verbosity
    fprintf('Apply integration over wavenumber domain ... '); 
end
% Apply quadrature weights.
u_FE = cellfun(@(x, y) {x * y}, u, num2cell(w));

% Sum up solutions.
u_FE = sum(cat(3, u_FE{:}), 3);

if verbosity
    fprintf('done.\n');
    fprintf('... DC 2.5D problem solved.\n \n');
end

% Obtain potentials at RX positions.
phi_FE = fe.I * u_FE;

%% Compare with analytic point source at top of homogeneous half-space.

if debugging
    % Get reference solution.
    u_ref = RefSol.getElectrodeAtHS(sig_anomaly, 2 * TX.val, TX.coo);
    phi_ref = arrayfun(@(x, y) u_ref.f(x, y), RX.coo(:, 1), RX.coo(:, 2));

    % Get asymptotic solution.
    phi_asy = zeros(size(RX.coo, 1), 1);
    for i=1:size(RX.coo, 1)
        r = abs(RX.coo(i,1) - TX.coo(1,1));
        phi_asy(i) = 2 / pi * sum(w .* besselk(0, k * r));
    end

    % Plot / compare solutions along the profile.
    x_plot = sqrt((RX.coo(:,1) - RX.coo(1,1)) .^2 + (RX.coo(:,2) - RX.coo(1,2)) .^2);
    figure()
    subplot(2, 1, 1)
        plot(x_plot, phi_ref, 'ok', x_plot, phi_asy, 'xb', x_plot, phi_FE, '-r');
        legend('\phi_{ref}', '\phi_{asy}', '\phi_{FE}');
        ylabel('potential');
    subplot(2, 1, 2)
        rel_err_FE = (1 - (phi_FE ./ phi_ref)) * 100;
        rel_err_asy = (1 - (phi_asy ./ phi_ref)) * 100;
        plot(x_plot, rel_err_FE, 'b', x_plot, rel_err_asy, 'r');
        legend('\phi_{ref} vs. \phi_{FE}', '\phi_{ref} vs. \phi_{asy}');
        ylabel('rel. error');
        xlabel('profile length');
end