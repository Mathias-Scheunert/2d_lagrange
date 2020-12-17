function mesh = NcreateNCVerticalDikeMesh(ele, x0, width,ref)
% Creates Gmsh input (*.geo) and output (.msh) file for a noncontigious
% vertical dike in a half space
% with threshold fields around the electrodes (y = 0).
%
% SYNTAX
%   [] = NcreateNCVerticalDikeMesh(ele, h, x0)
%
%     ___________________________
%     |       |       |          |
%     |       |       |          |
%     |       ---------          |
%     |                          |
%     |__________________________|
%
% INPUT PARAMETER
%   ele   ... Vector of electrode positions.
%   x0    ... Scalar, denoting left boundary of dike.
%   width ... Scalar, denoting width of dike.
%   ref   ... Scalar, denoting number of uniform grid refinements.
%
% OUTPUT PARAMETER
%   mesh ... Struct, containing vertex2coordinates (vertices),
%            simplex2vertex (cell2vtx), boundary egde2coordinates
%            (bnd_edge2vtx), boundary egde 2 model domain boundary
%            (bnd_edge_...) and simplex parameter domains
%            (parameter_domains).
%% Check input.

    if nargin < 4
        ref = 0;
    else
        assert(isscalar(ref));
    end

%% Check Gmsh version.

    gmsh_path = dir('**/gmsh.exe');
    if isempty(gmsh_path)
        error('Gmsh executable could not be found in path.');
    end
    %% Create Gmsh input file.

    gmsh_file = 'vert_dike';
    createGmshInput(ele, width, x0);

%% Run Gmsh.

    system([gmsh_path.folder, '/gmsh -2 ', ...
            gmsh_file, '.geo -v 0 -format msh2']);

    %% Import mesh information from .msh file.
    mesh = Mesh.loadGmsh([gmsh_file, '.msh'], ...
                         'ref', ref);
    % Override type.
    mesh.type = 'gmsh_load';

    %% Clean up.

    %delete([gmsh_file, '.geo'], [gmsh_file, '.msh']);

end

function createGmshInput(point_ele, h, x0)
    %  Note that the 2D grid will be referred to gmesh x and z coordinate:
    %  0,0 ------>
    %      |      x = gmesh-x
    %      |
    %      |
    %      v
    %       y = gmesh-z

    %% Domain boundary points.

    % If electrode positions are known, add an adapted offset.
    if ~isempty(point_ele)
        % Get specific extentions and adapt
        coo_offset = round(abs(max(point_ele(:,1)) - min(point_ele(:,1))));
        bnd = [-600*coo_offset, ...
               600*coo_offset, ...
               0, ...
               600*coo_offset];
        dike_depth =  600*coo_offset/4;
    end

    %middle = floor(min_x + coo_offset/2);% for potential secondary distance field

    % Transform boundary coords to boundary point coords.
    point_bnd = [bnd(1:2).', bnd(3)+zeros(2,1);
                 bnd(1:2).', bnd(4)+zeros(2,1)];
	n_point_bnd = size(point_bnd, 1);
    h_bnd = ['corner_mf'; 'corner_mf';'corner_mf';'corner_mf'];
    %% Dike points.

    point_dike = [x0, bnd(3);
                  x0 + h, bnd(3);
                  x0, dike_depth;
                  x0 + h, dike_depth];

    %% Electrode points.

    % Remove dike points which coincide with electrode positions.
    ele_at_dike = all(ismember(point_dike, point_ele).').';
    point_dike(ele_at_dike, :) = [];
    n_point_dike = size(point_dike, 1);
    n_point_ele = size(point_ele, 1);
    h_point_dike_ele = zeros((n_point_dike+n_point_ele), 1);
    h_point_dike_ele = string(h_point_dike_ele);

    %% Create point list input,

    % Summarize points.
    point_list = [point_bnd; ...
                  point_dike; ...
                  point_ele];
    n_point = n_point_bnd + n_point_dike + n_point_ele;
    point_h = [h_bnd; h_point_dike_ele];

    % Sort ascendingly.
    % (Note: only when sorted, line_list definition works!)
    [point_list, sort_idx] = sortrows(point_list);
    point_h = point_h(sort_idx);

    %  Create Gmsh [n x 5] matrix point-input structure.
    point_input = zeros(n_point, 4);
    for ii = 1:n_point
       point_input(ii, 1) = ii ;               % point number
       point_input(ii, 2) = point_list(ii, 1); % Gmsh-x cooridnate
       point_input(ii, 3) = 0;                 % Gmsh-y cooridnate
       point_input(ii, 4) = point_list(ii, 2); % Gmsh-z cooridnate
    end
    ele_search = sort_idx>= (length(sort_idx)+1-n_point_ele);
    ele_idx = point_input(ele_search,1);
    %% Handle geometric lines / line loops.

    % Get point indices in correct order (domain 1)
    pt_top_left = find(point_list(:, 1) <= x0 & point_list(:, 2) == 0);
    pt_top_right = find(point_list(:, 1) >= (x0 + h) & point_list(:, 2) == 0);
    pt_bot = flipud(find(point_list(:, 2) == bnd(4)));
    pt_bot_dike = find(point_list(:, 2) == dike_depth);

    idx_domain_1=[pt_top_left; pt_bot_dike; pt_top_right; pt_bot];
    line_1 = [idx_domain_1, [idx_domain_1(2:end); idx_domain_1(1)]];
    n_line_1 = size(line_1, 1);

    % Get point indices in correct order (dike).
    pt_top_mid = find((point_list(:, 1) >= x0 & point_list(:, 1) <= (x0 + h)) &...
                      point_list(:, 2) == 0);
    pt_bot_mid = flipud(find((point_list(:, 1) >= x0 & point_list(:, 1) <= (x0 + h)) &...
                              point_list(:, 2) == dike_depth)...
                       );
    idx_domain_2 = [pt_top_mid; pt_bot_mid];
    % Create line definition (domain 2).
    line_2 = [idx_domain_2, [idx_domain_2(2:end); idx_domain_2(1)]];
    n_line_2 = size(line_2, 1);

    % Create line input ([n x 3] matrix) line-structure.
    line_list = [line_1; line_2];
    n_line = n_line_1 + n_line_2;
    line_input = [(1:n_line).', line_list];

    %% Extract boundary line indices. (??

    pt_top = find(point_list(:, 2) == 0);
    pt_bot = find(point_list(:, 2) == bnd(4));
    pt_left = find(point_list(:, 1) == bnd(1));
    pt_right = find(point_list(:, 1) == bnd(2));

    idx_line_top = find(all(ismember(line_list, pt_top).'));
    idx_line_bot = find(all(ismember(line_list, pt_bot).'));
    idx_line_left = find(all(ismember(line_list, pt_left).'));
    idx_line_right = find(all(ismember(line_list, pt_right).'));

    %% Summarize output.

    fileID = fopen('vert_dike_params.geo', 'w');  %can easily be changed by other scripts
    fprintf(fileID, 'max_off = %f;\n',coo_offset);
    fprintf(fileID, 'dom = 600*max_off;\n');
    fprintf(fileID, 'dom_wdt = 2*dom;\n');
    fprintf(fileID, 'lc_min = 0.05;\n');
    fprintf(fileID, 'lc_max = dom_wdt/2;\n');
    fprintf(fileID, 'corner_mf = dom_wdt/2;\n');
    fprintf(fileID, 'DistMax = dom_wdt/0.5;\n');
    fprintf(fileID, 'DistMin = 0.8;\n');
    fclose(fileID);

    name = 'vert_dike.geo';
    fileID = fopen(name, 'w');
    fprintf(fileID,'// definitions \n');
    fprintf(fileID,'Include "vert_dike_params.geo"; \n');
    fprintf(fileID,'Mesh.Algorithm=1; \n');

   %% Handle geometric entities.
   % Points
   for ii = 1:n_point
        fprintf(fileID, 'Point(%d) = {%d, %d, %d, %s };\n', ...
            point_input(ii,:),point_h(ii,:));
   end

   % Add lines.
    fprintf(fileID, '\n');
    for ii = 1:n_line
        fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
            line_input(ii,:));
    end

    % Add line loops (polygonal chains).
    fprintf(fileID, '\n');
    fprintf(fileID, ...
            ['Line Loop(%d) = {', repmat('%d, ', 1, n_line_1-1), '%d};\n'], ...
            1, 1:n_line_1);
    fprintf(fileID, ...
            ['Line Loop(%d) = {', repmat('%d, ', 1, n_line_2-1), '%d};\n'], ...
            2, (n_line_1+1):(n_line_1+n_line_2));

    % Add geomentric surface.
    for ii = 1:2
        fprintf(fileID, ...
            'Plane Surface(%d) = {%d};\n', ...
            ii, ii);
    end

    %% Handle physical entities.
    fprintf(fileID, '\n');
    bnd_name = {'ymin', 'ymax', 'xmin', 'xmax'};
    idx_bnd = {idx_line_top, idx_line_bot, idx_line_left, idx_line_right}.';
    for kk = 1:length(bnd_name)
        fprintf(fileID, ...
            ['Physical Line("%s") = {', ...
                repmat('%d, ', 1, length(idx_bnd{kk})-1), '%d};\n'], ...
            bnd_name{kk}, idx_bnd{kk});
    end

    % Add physical surfaces.
    fprintf(fileID, '\n');
    dom_name = {'res_1', 'res_2'};
    for ii = 1:2
        fprintf(fileID, ...
            'Physical Surface("%s") = {%d};\n', ...
            dom_name{ii}, ii);
    end

    % Add Fields
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
            2*(1:(length(ele_idx))));
    fprintf(fileID,'Background Field = {%d};\n', 2*length(ele_idx)+1);

    % Remove multiples.
    fprintf(fileID, '\n');
    fprintf(fileID, 'Coherence;');

    fclose(fileID);
end



