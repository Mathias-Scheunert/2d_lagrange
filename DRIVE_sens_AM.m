% Analytische und numerische Sensitivitäten (bzgl. rho!)
%
% Sensitivitäten sind damit unabhängig vom Hintergrundwiderstand.

clearvars;
verbosity = pick(1, false, true);

%% Skript params

rho_val = 3;    % Hintergrundwid.
I = 1;          % Quellstärke
ngl = 6;        % Quad.ord. für Friedel-Ansatz (max 6)
refinement = 1; % Homogene Gitterverfeinerungen
FE_order = pick(2, 1, 2);

% Elektrodenpositionen (x, y=z=0) der 4-Pkt.-Anordnung.
ABMN_ele = [0 22.5 7.5 15];

%%

% Define type of numerical integration approach.
FT_type = 'Boerner';

ABMN = [ABMN_ele(1),ABMN_ele(3)]; %AM
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
                 'ref', fwd_params.ref, 'verbosity', verbosity);
mesh.cell2surf = Mesh.getCellArea(mesh);
fprintf('.');

%% Set up parameter vector.

param = Param.initParam(mesh, param_info);

%% Set up.

TX = struct();
TX.val = I;
TX.type = 'point_exact';
fprintf('.');
TX.coo = [ABMN(1), 0];
RX.coo = [ABMN(2), 0];

% Summarize parameter.
fwd_params.TX = TX;
fwd_params.RX = RX;

%% Assemble 2.5D DC problem mit sol.TA

[fe, sol, FT_info, bnd_reduc] = App_DC.NassembleDC25D_s(mesh, param, fwd_params, verbosity);

%% Solving FT_info.n 2D fwp.

u_2D = cell(FT_info.n, 1);
u_2D_reduc = cell(FT_info.n, 1);
tensor_product = cell(FT_info.n, 1);
for jj = 1:(FT_info.n)
        u_2D{jj} = FeL.solveFwd(sol{jj}, fe);
        u_2D_reduc{jj} = u_2D{jj}(bnd_reduc);
        tensor_product{jj} = full(ttv(sol{jj}.TA, u_2D_reduc{jj},2));
end
fprintf('.');

%% Computing J by solving fwp with Q

Q = fe.I(bnd_reduc);
J = zeros(size(Q, 1), size(mesh.cell2vtx, 1));
for ii = 1:(FT_info.n)
    QAinv = (sol{ii}.A \ Q')';
    % Applying quadrature approach.
    J_temp = -(1/pi)*full(QAinv *  tensor_product{ii} * FT_info.w(ii));
    J = J + J_temp; % "Numerische" Sensitivität bzgl. sigma
end
fprintf('.');

%% Visualize Mesh and TX/RX positions.

% AM
J_FE = J./mesh.cell2surf.'; % "Physikalische" Sensitivität
J_FE = -J_FE.*param'.^2;      % bzgl. rho

% Plot
ylim_max = ((max(ABMN_ele)+2) - (min(ABMN_ele)-2)) / 2;
ylim_min = -0.5;
colrang = Plot.plotSens(mesh, full(J_FE), 1, ABMN_ele, ylim_min, ylim_max);
title('numeric');
drawnow;

%% Analytic Sens

ABMN_ele_ana = [ABMN_ele(1), ABMN_ele(3)];

% Compute sensitivity
J_ana = Sens.sensitivitaet_msh(mesh, ABMN_ele_ana, ngl)'; % bzgl. rho
%FIXME: Warum ein Offset 2?
fprintf('.');

% Plot
% Plot.plotSens(mesh, J_ana, 2, ABMN_ele, ylim_min, ylim_max, colrang);
Plot.plotSens(mesh, J_ana, 2, ABMN_ele, ylim_min, ylim_max);
title('analytic');
drawnow;

%% Comparison

J_err = (J_FE - J_ana);
% Plot
Plot.plotSens(mesh, J_err, 3, ABMN_ele, ylim_min, ylim_max);
title('abs. error');
drawnow;

fprintf(sprintf('\nJ_err ... min: %f, max: %f \n', min(J_err), max(J_err)));

%% Taylortest für J

% Initialisieren.
% rng(0815);
h = logspace(2, -8, 11);
Q_full = fe.I;
u_FE = App_DC.solveDC25D(fe, sol, FT_info);
u_FE = Q_full*u_FE;

% J wählen.
J_tt_num = J; % numerische J
J_tt_ana = -J_ana.*mesh.cell2surf.'./param'.^2; % analytische J -> numerische J

% Referenzlösung für homog. HR
u_hr = @(sigma) (I)./(sigma*2*pi*(ABMN(2) - ABMN(1)));

% no param trafo
% ds2dm = 1;
% m2c = @(c) c;
% c2m = @(m) m;

% m = log(param) - trafo
ds2dm = diag(param);
m2c = @(c) log(c);
c2m = @(m) exp(m);

% arbitrary perturbation
% delta_m = randn(size(param))*c2m(param);

% uniform perturbation
delta_m = (randn*c2m(param(1)) + zeros(size(param)));

% delta_m = delta_m / norm(delta_m);

% Lösung ausrechnen
if ~exist('u_FE_tt', 'var')
    u_FE_tt = zeros(length(h), 1);
    for ii = 1:length(h)
        if ii > 1
          for j = 0:log10(ii-1)
              fprintf('\b');
          end
        end
        fprintf('%d', ii);
        param_tt = c2m(param) + h(ii)*delta_m;
        [~, sol_tt, ~] = App_DC.assembleDC25D(mesh, m2c(param_tt), fwd_params);
        u_FE_tt(ii) = Q_full*App_DC.solveDC25D(fe, sol_tt, FT_info);
    end
    fprintf('\n');
end

% Normen ausrechnen.
[err_tay_num1, err_tay_num2, err_tay_ana1, err_tay_ana2] = deal(zeros(length(h), 1));
for ii = 1:length(h)
    %num
    % || u(m_0 + h*dm) - u(m_0)||_2
    err_tay_num1(ii) = norm(u_FE_tt(ii)-u_FE, 2);
    % || u(m_0 + h*dm) - [u(m_0) + h*J(m_0)*dm]||_2
    err_tay_num2(ii) = norm(u_FE_tt(ii)-(u_FE + J_tt_num*ds2dm*h(ii)*delta_m), 2);
    %semi-ana
    err_tay_ana1(ii) = norm(u_FE_tt(ii)-(u_FE + J_tt_ana*ds2dm*h(ii)*delta_m), 2);
    %ana
    err_tay_ana2(ii) = norm(u_hr(c2m(param(1))+h(ii)*delta_m(1)) - ...
                           (u_hr( c2m(param(1))) + J_tt_ana*ds2dm*h(ii)*delta_m), 2);
    fprintf('.');
end
fprintf('\n');

% Konvergenzplot
figure(10);
loglog(h, err_tay_num1/err_tay_num1(1), 'x-m', ...
       h, err_tay_num2/err_tay_num2(1), 'd-m', ...
       h, err_tay_ana1/err_tay_ana1(1), '+-r', ...
       h, err_tay_ana2/err_tay_ana2(1), 'o-r');
hold on
    loglog(h, h/h(1), '--k');
    loglog(h, h.^2/h(1).^2, '.k');
hold off
title('Taylor @ FE');
ylim([1e-17, 1e2]);
xlim([min(h), max(h)]);
legend('e_0^{num}(h)', 'e_1^{num}(h)', 'e_1^{ana-num}(h)', 'e_1^{ana}(h)', 'O(h)', 'O(h^2)', ...
       'Location', 'SouthWest');
set(gca, 'XDir', 'reverse');

% simple_taylor(param(1), delta_m(1), h);
