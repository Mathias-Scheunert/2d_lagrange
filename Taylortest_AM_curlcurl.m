% Analytische und numerische Sensitivitäten (bzgl. rho!)
%
% Sensitivitäten sind damit unabhängig vom Hintergrundwiderstand.
% Standard: 5m spacing
close all
% clearvars;
verbosity = pick(1, false, true);
% Define type of numerical integration approach.
FT_type = 'Boerner';

%% Skript params

rho_val = 10;  % Background resistivity for ana sensitivity, superfluous
I = 1;          % Quellstärke
ngl = 6;        % Quad.ord. für Friedel-Ansatz
refinement = 0; % Homogene Gitterverfeinerungen
FE_order = pick(2, 1, 2);

model = pick(1,'homo','dike','2layers','3layers','box','premade');

% Source and receiver locations at earth's surface.
AM_ele = [1 8];

TX = struct();
TX.coo = [AM_ele(1),0];
TX.val = I;
TX.type = 'point_exact';
RX = struct();
RX.coo = [AM_ele(2),0];

fprintf('.');


%% Set up boundary conditions.
% Note: ymin denotes earth's surface.
bnd.type = {'dirichlet', 'neumann'};
%         ymin ymax  xmin xmax
bnd.val = {{[];   0;  0;   0}, ...   % 1 for Dirichlet
           {0;  [];   [];  []}}; ... % 2 for Neumann
bnd.name = {'ymin', 'ymax', 'xmin', 'xmax'};
bnd.quad_ord = 1;

% Summarize parameter.
fwd_params = struct();
fwd_params.TX = TX;
fwd_params.RX = RX;
fwd_params.bnd = bnd;
fwd_params.FT_type = FT_type;
fwd_params.FE_order = FE_order;
fwd_params.ref = refinement;
fprintf('.');

%% Set up mesh
switch model
    case 'homo'
        param_info.val = 100;
        param_info.name = {'res_1'};
        mesh = Mesh.initMesh('gmsh_homo', ...
                     'TX', TX.coo, ...
                     'RX', RX.coo, ...
                     'ref', fwd_params.ref);

    case 'dike'
        param_info.val = [1/1000, 1/100]; % res halfspace, res dike
        param_info.name = {'res_1','res_2'};
        dike_info = [2.5, 2]; %start, width
        %dike_info = [-30.5, 30];

        mesh = Mesh.initMesh('gmsh_dike', ...
                     'TX', [TX.coo; RX.coo], ...
                     'dike', dike_info, 'ref', fwd_params.ref);
    case '2layers'
        param_info.val = [1/1000, 1/100];
        param_info.name = {'res_1','res_2'};
        thickness = 3;
        mesh = Mesh.Ncreate_two_layer(sortrows([TX.coo; RX.coo]), thickness ,verbosity, fwd_params.ref);

    case '3layers'
        param_info.val = [1/100, 1/1000,1/100]; % res halfspace, res dike
        param_info.name = {'res_1','res_2','res_3'};
        layers_y = [3,3.5];
        mesh = Mesh.create3Layers(sortrows([TX.coo; RX.coo]),layers_y, verbosity, fwd_params.ref);

    case 'box'
        box = [4,7,2,4]; % x1,x2,y1,y2
        write_ele_and_box(sort(AM_ele),box); %writes parameter file for box_3/box_4.geo
        if eq(loop,2)
            file_name = 'box_3ele';
        else
            file_name = 'box_4ele';
        end
        gmsh_path = dir('**/gmsh.exe');
        if isempty(gmsh_path)
            error('Gmsh executable could not be found in path.');
        end
        system([gmsh_path.folder,'/gmsh -2 ', file_name, '.geo -v 0 -format msh2']);
        mesh_type = 'gmsh_load';
        mesh = Mesh.initMesh(mesh_type, 'name', [file_name, '.msh'], ...
                 'ref', fwd_params.ref, 'verbosity', verbosity);
        param_info.val = [1/1000, 1/100];
        param_info.name = {'res_1','res_2'};
    case 'premade' %for premeshed files
        param_info.val = [1/100, 1/1000]; % adjust as needed
        param_info.name = {'res_1','res_2'};
        file_name = 'dummy';
        system(['/gmsh -2 ', file_name, '.geo -v 0 -format msh2']);
        mesh_type = 'gmsh_load';
        mesh = Mesh.initMesh(mesh_type, 'name', [file_name, '.msh'], ...
                 'ref', fwd_params.ref, 'verbosity', verbosity);
    otherwise
%         error('invalid model')
end

fprintf('.');

%% Set up parameter vector.

param = Param.initParam(mesh, param_info);
for it = 1: length(param) % Cell surfaces
    mesh.cell2surf(it) = abs((mesh.cell2cord{it}(1,1).*(mesh.cell2cord{it}(2,2)-mesh.cell2cord{it}(3,2))+mesh.cell2cord{it}(2,1).*(mesh.cell2cord{it}(3,2)-mesh.cell2cord{it}(1,2))+mesh.cell2cord{it}(3,1).*(mesh.cell2cord{it}(1,2)-mesh.cell2cord{it}(2,2)))/2);
end


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
J_FE = J; % "Numerische" Sensitivität bzgl. sigma
fprintf('.');
clear TX.coo RX.coo fwd_params.TX fwd_params.RX

J_FE = J_FE./mesh.cell2surf; %sign off k necessary if BA AB are flipped
J_FE = -J_FE.*param'.^2; % "Physikalische" Sensitivität bzgl. rho

%% Analytic Sens, only homogeneous

% Compute sensitivity
J_ana = 2*Sens.sensitivitaet_msh(mesh, [AM_ele(1), AM_ele(2)], ngl)'; % bzgl. rho
fprintf('.');

%% Taylortest
close all
clear u_FE_tt
% Initialize
rng(0815);
delta_m = randn(size(param));
%rand intervall a b
% a=-1;
% b=1;
% delta_m = a + (b-a).*rand(size(param));
%delta_m = (delta_m / norm(delta_m))*4;
h = logspace(1, -6, 10);
%h = logspace(1, -2, 41); %weiss
Q = fe.I;
Q_reduc = full(fe.I);
Q_reduc = Q_reduc(bnd_reduc);
u_FE = Q*App_DC.solveDC25D(fe, sol, FT_info);

% Compute disturbed solutions/potentials
if ~exist('u_FE_tt', 'var')
    u_FE_tt = zeros(length(u_FE), length(h));
    for ii = 1:length(h)
        param_tt = param + h(ii)*delta_m;
        param_first(ii) = param_tt(1); %#ok<SAGROW>
        [~, sol_tt, ~] = App_DC.assembleDC25D(mesh, param_tt, fwd_params);
        u_FE_tt(:, ii) = Q*App_DC.solveDC25D(fe, sol_tt, FT_info);
        fprintf('.');
    end
    fprintf('\n');
end

% choose J
J_tt_num = J; % numeric J
J_tt_ana = -J_ana.*mesh.cell2surf./param'.^2; % analytische J -> numerische J

% Compute norm
[err_tay_num1, err_tay_num2, err_tay_ana] = deal(zeros(length(h), 1));
for ii = 1:length(h)
    err_tay_num1(ii) = abs(u_FE_tt(:, ii)-u_FE);
    err_tay_num2(ii) = abs(u_FE_tt(:, ii)-(u_FE + J_tt_num*h(ii)*delta_m));
    err_tay_ana(ii) = abs(u_FE_tt(:, ii)-(u_FE + J_tt_ana*h(ii)*delta_m));

%     err_tay_num1(ii) = norm(u_FE_tt(:, ii)-u_FE, 2);
%     err_tay_num2(ii) = norm(u_FE_tt(:, ii)-(u_FE + J_tt_num*h(ii)*delta_m), 2);
%     err_tay_ana(ii) = norm(u_FE_tt(:, ii)-(u_FE + J_tt_ana*h(ii)*delta_m), 2);

    fprintf('.');
end
fprintf('\n');

figure(10);
loglog(h, err_tay_num1/err_tay_num1(1), 'x-m', ...
       h, err_tay_num2/err_tay_num2(1), 'd-g', ...
       h, err_tay_ana/err_tay_ana(1), 'o-k', ...
       h, h/h(1), '--b', ...
       h, h.^2/h(1).^2, '--r');
hold off
ylim([1e-17, 1e2]);
xlim([min(h), max(h)]);
legend('err_{0}', 'err_{num}', 'err_{ana}', 'h', 'h^2', ...
       'Location', 'SouthWest');
set(gca, 'XDir', 'reverse');
