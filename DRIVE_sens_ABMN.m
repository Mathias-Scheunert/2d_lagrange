% Analytische und numerische Sensitivitäten (bzgl. rho!)
%
% Sensitivitäten sind damit unabhängig vom Hintergrundwiderstand.

clearvars;
verbosity = pick(1, false, true);

%% Skript params

rho_val = 1;    % Hintergrundwid.
I = 1;          % Quellstärke
ngl = 6;        % Quad.ord. für Friedel-Ansatz
refinement = 1; % Homogene Gitterverfeinerungen
FE_order = pick(2, 1, 2);

% Elektrodenpositionen (x, y=z=0) der 4-Pkt.-Anordnung.
% ABMN_ele = [-1.5 1.5 -0.5 0.5];
ABMN_ele = [0 22.5 7.5 15];
% ABMN_ele = [0 99 33 66];

%%

% Define type of numerical integration approach.
FT_type = 'Boerner';

ABMN = [ABMN_ele(1),ABMN_ele(3); %AM
        ABMN_ele(1),ABMN_ele(4); %AN
        ABMN_ele(2),ABMN_ele(3); %BM
        ABMN_ele(2), ABMN_ele(4)]; %BN
loop = 4; %for for loop

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

TX = struct();
TX.val = I;
TX.type = 'point_exact';
fprintf('.');

%% Loop over electrodes
JJ = cell(loop,1);
for count = 1 : loop %4 pole pole  %2 bei pdp

    %% Set up
    TX.coo = [ABMN(count,1), 0];
    RX.coo = [ABMN(count,2), 0];

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
            tensor_product{jj} = full(ttv(sol{jj}.TA,u_2D_reduc{jj},2));
    end
    fprintf('.');

    %% Computing J by solving fwp with Q

    Q = fe.I(bnd_reduc);
    J = zeros(size(Q, 1), size(mesh.cell2vtx, 1));
    for ii = 1:(FT_info.n)
        QAinv = (sol{ii}.A \ Q')';
        J_temp = -1/pi * full(QAinv *  tensor_product{ii}*FT_info.w(ii));
        J = J + J_temp;
    end

    % Normalize to cell size.
    JJ{count} = J; % "Numerische" Sensitivität bzgl. sigma
    fprintf('.');
end

%% Visualize Mesh and TX/RX positions.

% (AM - AN) - (BM - BN)
J_FE = (JJ{1}-JJ{2})-(JJ{3}-JJ{4});
J_FE = J_FE./mesh.cell2surf.';
J_FE = -J_FE.*param'.^2; % "Physikalische" Sensitivität bzgl. rho

% Plot
ylim_max = 10;
ylim_min = -1;
colrang = Plot.plotSens(mesh, full(J_FE), 1, ABMN_ele, ylim_min, ylim_max);
title('numeric');
drawnow;

%% Analytic Sens

ABMN_ele = [ABMN_ele(1), ABMN_ele(3), ABMN_ele(4), ABMN_ele(2)];

% Compute sensitivity
J_ana = Sens.sensitivitaet_msh(mesh, ABMN_ele, ngl)'; % bzgl. rho
fprintf('.');

% Plot
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
