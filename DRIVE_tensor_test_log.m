% Analytische und numerische Sensitivitäten (bzgl. rho!)
%
% Sensitivitäten sind damit unabhängig vom Hintergrundwiderstand.

clearvars;

%% Skript params.

rho_val = 1;    % Hintergrundwid.
I = 1;          % Quellstärke
ngl = 6;        % Quad.ord. für Friedel-Ansatz (max 6)
refinement = 0; % Homogene Gitterverfeinerungen
FE_order = pick(2, 1, 2);
FT_type = 'Boerner';

% Elektrodenpositionen (x, y=z=0) der 4-Pkt.-Anordnung.
ABMN_ele = [0 22.5 7.5 15];
fprintf('.');

%% Set up boundary conditions.

% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax
bnd.val = {{[];   0;  0;   0}, ...   % 1 for Dirichlet
           {0;  [];   [];  []}}; ... % 2 for Neumann
bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};
bnd.quad_ord = 1;

% Define background conductivity
param_info = struct();
param_info.val = 1/rho_val;
param_info.name = {'res_1'};

fwd_params = struct();
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;
clear('TX', 'RX', 'bnd', 'FT_type', 'FE_order', ...
      'refinement', 'topo', 'x', 'y');
fprintf('.');

%% Set up mesh.

file_name = '+Sens/default';
Sens.write_pt(ABMN_ele);
system(['gmsh -2 ', file_name, '.geo -v 0 -format msh2']);
mesh_type = 'gmsh_load';
mesh = Mesh.initMesh(mesh_type, 'name', [file_name, '.msh'], ...
                                'ref', fwd_params.ref);
mesh.cell2surf = Mesh.getCellArea(mesh);
fprintf('.');

%% Set up parameter vector.

param = Param.initParam(mesh, param_info);

%% Set up.

TX = struct();
TX.val = I;
TX.type = 'point_exact';
fprintf('.');
TX.coo = [ABMN_ele(1), 0];
RX.coo = [ABMN_ele(3), 0];

% Summarize parameter.
fwd_params.TX = TX;
fwd_params.RX = RX;

%% Assemble 2.5D DC problem mit sol.TA

s2m = @(x) log(x);
m2s = @(x) exp(x);
m = s2m(param);

[fe, sol, FT_info, bnd_reduc] = App_DC.NassembleDC25D_s(mesh, m2s(m), fwd_params);
fprintf('.');
fprintf('\n');

%% Taylortest für Tensor.

% Initialisieren.
rng(0815);
d_m = randn(size(m));
d_m = d_m / norm(d_m);
h = logspace(0, -8, 51);
n_h = length(h);

% Über alle h iterieren und A(h) sowie TA(h) assemblieren.
sol_tt = cell(length(sol), n_h);
for ii = 1:n_h
    % Neues Problem assemblieren.
    m_tt = m + h(ii)*d_m;
    [~, tmp, ~, ~] = App_DC.NassembleDC25D_s(mesh, m2s(m_tt), fwd_params);
    sol_tt(:, ii) = tmp;
    fprintf('.');
end
fprintf('\n');

%% Normen berechnen und darstellen.

% Über alle Wellenzahlen iterieren.
for kk = 1:FT_info.n
    [err_tt_1, err_tt_2, err_tt_3] = deal(zeros(n_h, 1));
    for ii = 1:n_h
        % Normen ausrechnen.
        % || A(m_0 + h*dm) - A(m_0)||_fro
        err_tt_1(ii) = norm(sol_tt{kk, ii}.A - sol{kk}.A, 'fro');
        % || A(m_0 + h*dm) - [A(m_0) + h*(TA(m_0) x_3 dm)]||_fro
        err_tt_2(ii) = norm(sol_tt{kk, ii}.A - ...
                           (sol{kk}.A + ttv(sol{kk}.TA, h(ii)*diag(m2s(m))*d_m, 3)), ...
                            'fro');
        % Nur Tensor(en) verwenden.
        A_TA = ttv(sol{kk}.TA, m2s(m + h(ii)*d_m), 3);
        A_TA0 = ttv(sol{kk}.TA, m2s(m), 3);
        err_tt_3(ii) = norm(A_TA - A_TA0 - ttv(sol{kk}.TA, h(ii)*diag(m2s(m))*d_m, 3), ...
                            'fro');
    end

    % Konvergenzplot
    figure(10);
    loglog(h, err_tt_1/err_tt_1(1), 'x-m', ...
           h, err_tt_2/err_tt_2(1), 'o-m', ...
           h, err_tt_3/err_tt_3(1), 'd-r');
    hold on
        loglog(h, h/h(1), '--k');
        loglog(h, h.^2/h(1).^2, '.k');
    hold off
    title(sprintf('Taylor @ A(k = %f)', FT_info.k(kk)));
    ylim([1e-17, 1e2]);
    xlim([min(h), max(h)]);
    legend('e_0(h)', 'e_1(h)', 'e_{1 test}(h)', 'O(h)', 'O(h^2)', ...
           'Location', 'SouthWest');
    set(gca, 'XDir', 'reverse');
    pause(0.15)
end
