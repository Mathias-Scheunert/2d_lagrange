function mesh = Ncreate_two_layer(ele, thickness ,verbosity,ref)
    %Creates mesh with two resistivity domains

    %% Check Gmsh version.

    gmsh_path = dir('**/gmsh.exe');
    if isempty(gmsh_path)
        error('Gmsh executable could not be found in path.');
    end

    %% Create Gmsh input file.

    gmsh_file = 'two_layers';
    domains = {'res_1','res_2'};
    createGeoFile([gmsh_file, '.geo'], ele, domains, thickness, verbosity);

    %% Run Gmsh.

    % TODO: if it occurs, catch error (message) and stop.
    %     system([gmsh_path.folder, '/gmsh -2 ', ...
    %             gmsh_file, '.geo -v 0 -format msh2']);
     system([gmsh_path.folder, '/gmsh -2 ', ...
            gmsh_file, '.geo -v 0 -format msh2 ']);

    %% Import mesh information from .msh file.
    mesh_type = 'gmsh_load';
    mesh = Mesh.initMesh(mesh_type, 'name', [gmsh_file, '.msh'], ...
                 'ref', ref, 'verbosity', verbosity);
%     mesh = Mesh.loadGmsh([gmsh_file, '.msh'], ...
%         'verbosity', verbosity, ...
%         'ref', ref);
%
%     mesh.type = 'gmsh_load';
end

function createGeoFile(name, ele, dom, thickness , verbosity)
     if verbosity
       fprintf('Define Gmsh input ... ');
     end
    %% Set domain boundaries.

    TXRX = unique(ele, 'rows');
    if ~isempty(TXRX)

        % Get specific extentions.
        min_x = min(TXRX(:,1));
        max_x = max(TXRX(:,1));
        offset = max_x-min_x;

        % Adapt.
        bnd = [min_x - offset*600, ...
               max_x + offset*600, ...
               0, ...
               offset*600];
    end
% Transform boundary coords to boundary point coords.
    point_bnd = [bnd(1:2).', bnd(3)+zeros(2,1);
                 bnd(1:2).', bnd(4)+zeros(2,1)];
%% Set up domain basic geometry entities (points, lines).

    point_domain = zeros(6, 5);
    for i = 1:4
       point_domain(i, 1) = i ;               % point number
       point_domain(i, 2) = point_bnd(i,1);   % gmesh-x cooridnate.
       point_domain(i, 3) = 0;                % gmesh-y cooridnate.
       point_domain(i, 4) = point_bnd(i,2);   % gmesh-z cooridnate.
       point_domain(i, 5) = 0;
    end
    point_domain(3:4, 5) = offset*600/2;

    point_domain(5,:) = [5,point_bnd(1,1),0,thickness,0];
    point_domain(6,:) = [6,point_bnd(2,1),0,thickness,0];

    %% Electrode points.
    point_ele = zeros(size(ele,1),5);
    for i = 1:size(ele,1)
       point_ele(i, 1) = i + 6;      % point number
       point_ele(i, 2) = ele(i,1);   % gmesh-x cooridnate.
       point_ele(i, 3) = 0;          % gmesh-y cooridnate.
       point_ele(i, 4) = ele(i,2);   % gmesh-z cooridnate.
       point_ele(i, 5) = 0;
    end
    % Summarize points.
    point_list = [point_domain; ...
                  point_ele];
    ele_idx = 7:(size(ele,1)+6);
    %[point_list, sort_idx] = sortrows(point_list,2);

    %% Handle geometric lines / line loops.
    lines_bnd = [1,5; 5,3; 5,6; 3,4; 2,6; 6,4;]; %2left, middle, down, 2right
    lines_ele = zeros(size(ele,1)+1,2);
    lines_ele(1,:) = [1,7];
    for i = 1:(size(ele,1)-1)
       lines_ele(i+1,:) = [point_ele(i,1),point_ele(i+1,1)];
    end
    lines_ele(size(ele,1)+1,:)= [size(point_list,1),2];

     % Create line input ([n x 3] matrix) line-structure.
    line_list =[lines_bnd; lines_ele];
    line_list = [(1:size(line_list,1))',line_list];

    %Loops
    loop{1} = [1,3,-5,-flip(line_list(7:end,1))'];
    loop{2} = [2,4,-6,-3];

    dom_wdt = (abs(bnd(1))+abs(bnd(2)))/2;

    %%Write .geo
    fileID = fopen(name, 'w');
    fprintf(fileID,'// definitions \n');
	fprintf(fileID, 'dom_wdt = %f;\nlc_min = 0.05;\nlc_max = dom_wdt/2;\nDistMax = dom_wdt/0.5;\nDistMin = 0.8;\n Mesh.Algorithm=1; \n',dom_wdt);

    % Add points.
    fprintf(fileID, '\n');
    for ii = 1:size(point_list,1)
        fprintf(fileID, 'Point(%d) = {%d, %d, %d, %d};\n', ...
            point_list(ii,:));
    end

   % Add lines.
    fprintf(fileID, '\n');
    for ii = 1:size(line_list,1)
        fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
            line_list(ii,:));
    end

    % Add line loops (polygonal chains).
    fprintf(fileID, '\n');
    for ii = 1:2
        cur_n_line = size(loop{ii}, 2);
        fprintf(fileID, ...
            ['Line Loop(%d) = {', repmat('%d, ', 1, cur_n_line-1), '%d};\n'], ...
            ii, loop{ii});
    end
    fprintf(fileID, '\n');
    for ii = 1:2
        fprintf(fileID, ...
            'Plane Surface(%d) = {%d};\n', ...
            ii, ii);
    end
%% Handle physical entities.
    fprintf(fileID, '\n');
    % Add physical lines.
    fprintf(fileID, '\n');
    fprintf(fileID, ...
	'Physical Line("%s") = {%d,%d};\n','xmin', [1,2]);
    fprintf(fileID, ...
	'Physical Line("%s") = {%d,%d};\n','xmax', [5,6]);
    fprintf(fileID, ...
	'Physical Line("%s") = {%d};\n','ymax', 4);
    n_lines_top = size(line_list(7:end,:),1);
    fprintf(fileID, ...
    ['Physical Line("%s") = {', ...
        repmat('%d, ', 1, n_lines_top-1), '%d};\n'], ...
    'ymin', line_list(7:end,1));

    % Add physical surfaces.
    fprintf(fileID, '\n');
    for ii = 1:length(loop)
        fprintf(fileID, ...
            'Physical Surface("%s") = {%d};\n', ...
            dom{ii}, ii);
    end

%% Add Fields
    fprintf(fileID, '\n');
    for ii = 1:length(ele_idx)
        fprintf(fileID,'Field[%d] = Distance;\n', 2*ii-1);
        fprintf(fileID,'Field[%d].NodesList = {%d};\n', 2*ii-1,ele_idx(ii));
        fprintf(fileID,'Field[%d] = Threshold;\n', 2*ii);
        fprintf(fileID,'Field[%d].IField = %d;\n', 2*ii, 2*ii-1);
        fprintf(fileID,'Field[%d].LcMin = lc_min;\n', 2*ii);
        fprintf(fileID,'Field[%d].LcMax = lc_max;\n', 2*ii);
        fprintf(fileID,'Field[%d].DistMin = DistMin;\n', 2*ii);
        fprintf(fileID,'Field[%d].DistMax = DistMax;\n \n', 2*ii);
    end

    fprintf(fileID,'Field[%d] = Min;\n', 2*length(ele_idx)+1);
    fprintf(fileID,'Field[%d].FieldsList = {',2*length(ele_idx)+1);
    fprintf(fileID, [ ...
            repmat('%d, ', 1, length(ele_idx)-1), '%d};\n'], ...
            2*(1:length(ele_idx)));
    fprintf(fileID,'Background Field = {%d};\n', 2*length(ele_idx)+1);



    fclose(fileID);

end
