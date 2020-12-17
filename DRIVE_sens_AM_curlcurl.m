% Analytische und numerische Sensitivitäten (bzgl. rho!)
%
% Sensitivitäten sind damit unabhängig vom Hintergrundwiderstand.

clearvars;
verbosity = pick(1, false, true);

%% Skript params

rho_val = 1/1000;    % Hintergrundwid.
I = 1;          % Quellstärke
ngl = 6;        % Quad.ord. für Friedel-Ansatz (max 6)
refinement = 1; % Homogene Gitterverfeinerungen
FE_order = pick(2, 1, 2);

% Elektrodenpositionen (x, y=z=0) der 4-Pkt.-Anordnung.
ABMN_ele = [0 7.5];


x_min = min(ABMN_ele);
x_max = max(ABMN_ele);
pad = (x_max-x_min)/4;
x_lim = [x_min-pad, x_max+pad];
y_lim = [-2*pad, 1];

%%

% Define type of numerical integration approach.
FT_type = 'Boerner';

ABMN = [ABMN_ele(1),ABMN_ele(2)]; %AM
fprintf('.');

%% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax
bnd.val = {{[]; 0}, ...   % 1 for Dirichlet
           { 0; []}}; ... % 2 for Neumann
bnd.name = {'surface', 'subsurface'};
bnd.quad_ord = 1;

% Define background conductivity
param_info = struct();
param_info.val = 1/rho_val;
param_info.name = {'domain'};

fwd_params = struct();
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;
clear('TX', 'RX', 'bnd', 'FT_type', 'FE_order', ...
      'refinement', 'topo', 'x', 'y');
fprintf('.');

%% Set up mesh.

file_name = 'tmp_mesh_ABMN';
system(['gmsh -save ', file_name, '.geo -v 0 -format msh2']);
mesh_type = 'gmsh_load';
mesh_2d = Mesh.initMesh(mesh_type, 'name', [file_name, '.msh'], ...
                        'ref', fwd_params.ref, 'verbosity', verbosity);
mesh_2d.cell2surf = Mesh.getCellArea(mesh_2d);
fprintf('.');

%% Set up parameter vector.

param = Param.initParam(mesh_2d, param_info);

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

[fe_2d, sol_2d, FT_info, bnd_reduc] = App_DC.NassembleDC25D_s(mesh_2d, param, fwd_params, verbosity);

%% Solving FT_info.n 2D fwp.

u_2D = cell(FT_info.n, 1);
u_2D_reduc = cell(FT_info.n, 1);
tensor_product = cell(FT_info.n, 1);
for jj = 1:(FT_info.n)
        u_2D{jj} = FeL.solveFwd(sol_2d{jj}, fe_2d);
        u_2D_reduc{jj} = u_2D{jj}(bnd_reduc);
        tensor_product{jj} = full(ttv(sol_2d{jj}.TA, u_2D_reduc{jj},2));
end
fprintf('.');

%% Computing J by solving fwp with Q

Q = fe_2d.I(bnd_reduc);
J = zeros(size(Q, 1), size(mesh_2d.cell2vtx, 1));
for ii = 1:(FT_info.n)
    QAinv = (sol_2d{ii}.A \ Q')';
    % Applying quadrature approach.
    J_temp = -(1/pi)*full(QAinv *  tensor_product{ii} * FT_info.w(ii));
    J = J + J_temp; % "Numerische" Sensitivität bzgl. sigma
end
fprintf('.');

%% Visualize Mesh and TX/RX positions.

% AM
J_FE_2d = J./mesh_2d.cell2surf.'; % "Physikalische" Sensitivität
J_FE_2d = -J_FE_2d.*param'.^2;      % bzgl. rho

% Plot
ylim_max = 10;
ylim_min = 1;

colrang = Plot.plotSens(mesh_2d, full(J_FE_2d), 10, ABMN_ele, ylim_min, ylim_max);
xlim(x_lim);
ylim(y_lim);
title('numeric');
drawnow;

%% Analytic Sens

ABMN_ele_ana = [ABMN_ele(1), ABMN_ele(2)];

% Compute sensitivity
J_ana = 2*Sens.sensitivitaet_msh(mesh_2d, ABMN_ele_ana, ngl)'; % bzgl. rho
%FIXME: Warum ein Offset 2?
fprintf('.');

% Plot
% Plot.plotSens(mesh, J_ana, 2, ABMN_ele, ylim_min, ylim_max, colrang);
Plot.plotSens(mesh_2d, J_ana, 20, ABMN_ele, ylim_min, ylim_max);
xlim(x_lim);
ylim(y_lim);
title('analytic');
drawnow;

%% Comparison

J_err = (J_FE_2d - J_ana);
% Plot
Plot.plotSens(mesh_2d, J_err, 30, ABMN_ele, ylim_min, ylim_max);
xlim(x_lim);
ylim(y_lim);
title('abs. error');
drawnow;

fprintf(sprintf('\nJ_err ... min: %f, max: %f \n', min(J_err), max(J_err)));

%% Taylortest für J

% Initialisieren.
% rng(0815);
h = logspace(0, -8, 9);
Q_full = fe_2d.I;
u_FE = App_DC.solveDC25D(fe_2d, sol_2d, FT_info);
u_FE = Q_full*u_FE;

% J wählen.
J_ana2num = -J_ana.*mesh_2d.cell2surf.'./param'.^2; % analytische J -> numerische J

% Referenzlösung für homog. HR
u_hr = @(sigma) (I)./(sigma*2*pi*abs(ABMN(2) - ABMN(1)));

% no param trafo
dc2dm = 1;
c2m = @(c) c;
m2c = @(m) m;

% m = log(param) - trafo
% dc2dm = diag(param);
% c2m = @(c) log(c);
% m2c = @(m) exp(m);

m_0 = c2m(param);

% arbitrary perturbation
% delta_m = randn(size(param)).*m_0;

% uniform perturbation
delta_m = (randn*m_0(1) + zeros(size(param)));

% delta_m = delta_m / norm(delta_m);

% Lösung ausrechnen
u_FE_tt = zeros(length(h), 1);
for ii = 1:length(h)
    if ii > 1
      for j = 0:log10(ii-1)
          fprintf('\b');
      end
    end
    fprintf('%d', ii);
    m_tt = m_0 + h(ii)*delta_m;
    [~, sol_tt, ~] = App_DC.assembleDC25D(mesh_2d, m2c(m_tt), fwd_params);
    u_FE_tt(ii) = Q_full*App_DC.solveDC25D(fe_2d, sol_tt, FT_info);
end
fprintf('\n');


% Normen ausrechnen.
[err_tay_num0, err_tay_num1, err_tay_ana1] = deal(zeros(length(h), 1));
for ii = 1:length(h)
    %num
    % || u(m_0 + h*dm) - u(m_0)||_2
    err_tay_num0(ii) = norm(u_FE_tt(ii)-u_FE, 2);
    % || u(m_0 + h*dm) - [u(m_0) + h*J(m_0)*dm]||_2
    err_tay_num1(ii) = norm(u_FE_tt(ii)-(u_FE + J*dc2dm*h(ii)*delta_m), 2);
    %ana
    err_tay_ana1(ii) = norm(u_hr(m2c(m_0(1)+h(ii)*delta_m(1))) - ...
                           (u_hr(param(1)) + J_ana2num*dc2dm*h(ii)*delta_m), 2);
    fprintf('.');
end
fprintf('\n');

% Konvergenzplot
figure(200);
loglog(h, err_tay_num0/err_tay_num0(1), 'o-b', ...
       h, err_tay_num1/err_tay_num1(1), 'd-r', ...
       h, err_tay_ana1/err_tay_ana1(1), '+-r');
hold on
    loglog(h, h/h(1), '.b');
    loglog(h, h.^2/h(1).^2, '.r');
hold off
title('Taylor @ FE');
ylim([1e-17, 1e2]);
xlim([min(h), max(h)]);
legend('e_0^{num}(h)', 'e_1^{num}(h)', 'e_1^{ana}(h)', 'O(h)', 'O(h^2)', ...
       'Location', 'SouthWest');
set(gca, 'XDir', 'reverse');

% simple_taylor(param(1), delta_m(1), h);
