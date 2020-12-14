% Analytische und numerische Sensitivitäten (bzgl. rho!)
%
% Sensitivitäten sind damit unabhängig vom Hintergrundwiderstand.

clearvars;
verbosity = pick(1, false, true);

%% Skript params

rho_val = 1000;    % Hintergrundwid.
I = 1;          % Quellstärke
ngl = 6;        % Quad.ord. für Friedel-Ansatz
refinement = 0; % Homogene Gitterverfeinerungen
FE_order = pick(2, 1, 2);

% Elektrodenpositionen (x, y=z=0) der 4-Pkt.-Anordnung.
ABMN_ele = [0 22.5 7.5 15];

%%

% Define type of numerical integration approach.
FT_type = 'Boerner';

ABMN = [ABMN_ele(1),ABMN_ele(3);   %AM
        ABMN_ele(1),ABMN_ele(4);   %AN
        ABMN_ele(2),ABMN_ele(3);   %BM
        ABMN_ele(2), ABMN_ele(4)]; %BN
loop = 4; %for for loop

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
system(['gmsh  -save ', file_name, '.geo -v 0 -format msh2']);
mesh_type = 'gmsh_load';
mesh = Mesh.initMesh(mesh_type, 'name', [file_name, '.msh'], ...
                     'ref', fwd_params.ref, 'verbosity', verbosity);
mesh.cell2surf = Mesh.getCellArea(mesh);
fprintf('.');

%% Set up parameter vector.

param = Param.initParam(mesh, param_info);

TX = struct();
TX.val = I;
TX.type = 'point_exact';
fprintf('.');

%% Loop over electrodes
[JJ, DD] = deal(cell(loop,1));
for count = 1 : loop %4 pole pole  %2 bei pdp

    %% Set up
    TX.coo = [ABMN(count,1), 0];
    RX.coo = [ABMN(count,2), 0];

    % Summarize parameter.
    fwd_params.TX = TX;
    fwd_params.RX = RX;

    %% Assemble 2.5D DC problem mit sol.TA

    [fe, sol, FT_info, bnd_reduc] = App_DC.NassembleDC25D_s(mesh, param, fwd_params, verbosity);

    % Calculate J.
    Q = fe.I(bnd_reduc);
    J = zeros(size(Q, 1), size(mesh.cell2vtx, 1));
    for jj = 1:(FT_info.n)
        % Solving FT_info.n 2D fwp.
        u_2D = FeL.solveFwd(sol{jj}, fe);
        u_2D_reduc = u_2D(bnd_reduc);
        tensor_product = full(ttv(sol{jj}.TA,u_2D_reduc,2));

        % Computing J by solving fwp with Q
        QAinv = (sol{jj}.A \ Q')';
        J_tmp = -full(QAinv * tensor_product);
        J_tmp = 1/pi * J_tmp * FT_info.w(jj);
        J = J + J_tmp;
    end

    % Normalize to cell size.
    JJ{count} = J; % "Numerische" Sensitivität bzgl. sigma
    fprintf('.');

    %% Get observation.

    UU = App_DC.solveDC25D(fe, sol, FT_info);
    DD{count} = fe.I * UU;
end

%% Visualize Mesh and TX/RX positions.

% (AM - AN) - (BM - BN)
Neumann = @(x) (x{1}-x{2})-(x{3}-x{4});
D_FE = Neumann(DD);
J_FE_pre = Neumann(JJ);
J_FE = J_FE_pre;
J_FE = J_FE./mesh.cell2surf.';
J_FE = -J_FE.*param'.^2; % "Physikalische" Sensitivität bzgl. rho

% Plot
ylim_max = 1;
ylim_min = -10;
colrang = Plot.plotSens(mesh, full(J_FE), 10, ABMN_ele, ylim_min, ylim_max);
title('numeric');
drawnow;

%% Analytic Sens

AMNB_ele = [ABMN_ele(1), ABMN_ele(3), ABMN_ele(4), ABMN_ele(2)];

% Compute sensitivity
J_ana = 2*Sens.sensitivitaet_msh(mesh, AMNB_ele, ngl)'; % bzgl. rho
fprintf('.');

% Plot
Plot.plotSens(mesh, J_ana, 20, AMNB_ele, ylim_min, ylim_max);
title('analytic');
drawnow;

%% Comparison

J_err = (J_FE - J_ana);
% Plot
Plot.plotSens(mesh, J_err, 30, AMNB_ele, ylim_min, ylim_max);
title('abs. error');
drawnow;

fprintf(sprintf('\nJ_err ... min: %f, max: %f \n', min(J_err), max(J_err)));
