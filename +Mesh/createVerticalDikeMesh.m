function mesh = createVerticalDikeMesh(bnd, ele, x0, width, ref)
    % Creates Gmsh input (*.geo) and output (.msh) file for a VD in a HS.
    %
    % Constructs mesh for a half-space including a vertical dike and a
    % arbitrary number of TX and RX positions on top (i.e. at y = 0).
    %
    % SYNTAX
    %   [] = createVerticalDikeMesh(bnd, ele, h, x0[, gmsh_file])
    %
    % INPUT PARAMETER
    %   bnd   ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
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
    
    assert(all(ele(:,2) == 0) &&  bnd(3) == 0, ...
        'Only flat surface with electrodes at y = 0, supported.')
    assert(bnd(1) < x0 && bnd(2) > (x0 + width), ...
        'Expect dike to be located inside the given domain boundary.')
    if nargin < 5
        ref = 0;
    else
        assert(isscalar(ref));
    end
    
    %% Check Gmsh version.
    
    gmsh_path = dir('**/gmsh');
    if isempty(gmsh_path)
        error('Gmsh executable could not be found in path.');
    end
    
    %% Create Gmsh input file.

    gmsh_file = 'vert_dike';
    input_info = createGmshInput(bnd, ele, width, x0);
    writeGEO([gmsh_file, '.geo'], input_info)
    
    %% Run Gmsh.
    
    system([gmsh_path.folder, '/gmsh -2 ', ...
            gmsh_file, '.geo -v 0 -format msh2']);
        
    %% Import mesh information from .msh file.

    mesh = Mesh.loadGmsh([gmsh_file, '.msh'], ...
                         'ref', ref);
    
    % Override type.
    mesh.type = 'gmsh_load';

    %% Clean up.
    
    delete([gmsh_file, '.geo'], [gmsh_file, '.msh']);
end

function info = createGmshInput(bnd, point_ele, h, x0)
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
        % Get specific extentions.
        min_x = min(point_ele(:,1));
        max_x = max(point_ele(:,1));
        max_y = max(point_ele(:,2));
        
        % Adapt.
        coo_offset = round(abs(max_x - min_x) * 500);
        bnd = [floor(min_x - coo_offset), ...
               ceil(max_x + coo_offset), ...
               bnd(3), ...
               ceil(max_y + coo_offset)];
    end
    
    % Set maximal element size.
    h_bnd = round(5 * diff(bnd(1:2)) / 25, 1, 'significant'); % domain
    h_ele = 1;                                                % electrodes
    
    % Transform boundary coords to boundary point coords.
    point_bnd = [bnd(1:2).', bnd(3)+zeros(2,1);
                 bnd(1:2).', bnd(4)+zeros(2,1)];
	n_point_bnd = size(point_bnd, 1);

    %% Dike points.
    
    point_dike = [x0, bnd(3);
                  x0 + h, bnd(3);
                  x0, bnd(4);
                  x0 + h, bnd(4)];
    h_dike = [zeros(2, 1) + h_ele; zeros(2, 1) + h_bnd];
             
    %% Electrode points.

    % Remove dike points which coincide with electrode positions.
    ele_at_dike = all(ismember(point_dike, point_ele).').';
    point_dike(ele_at_dike, :) = [];
    h_dike(ele_at_dike, :) = [];
    n_point_dike = size(point_dike, 1);
    n_point_ele = size(point_ele, 1);
        
    %% Create point list input,
    
    % Summarize points.
    point_list = [point_bnd; ...
                  point_dike; ...
                  point_ele];
    point_h = [zeros(n_point_bnd, 1) + h_bnd; ...
               h_dike; ...
               zeros(n_point_ele, 1) + h_ele];
    n_point = n_point_bnd + n_point_dike + n_point_ele;
    
    % Sort ascendingly.
    % (Note: only when sorted, line_list definition works!)
    [point_list, sort_idx] = sortrows(point_list);
    point_h = point_h(sort_idx);
    
    %  Create Gmsh [n x 5] matrix point-input structure.
    point_input = zeros(n_point, 5);
    for ii = 1:n_point
       point_input(ii, 1) = ii ;               % point number
       point_input(ii, 2) = point_list(ii, 1); % Gmsh-x cooridnate
       point_input(ii, 3) = 0;                 % Gmsh-y cooridnate
       point_input(ii, 4) = point_list(ii, 2); % Gmsh-z cooridnate
       point_input(ii, 5) = point_h(ii);       % element size
    end    
    
    %% Handle geometric lines / line loops.
    
    % Get point indices in correct order (domain 1).
    pt_top_left = find(point_list(:, 1) <= x0 & point_list(:, 2) == 0);
    pt_bot_left = flipud(find(point_list(:, 1) <= x0 & ...
                              point_list(:, 2) == bnd(4)) ...
                         );
    idx_domain_1 = [pt_top_left; pt_bot_left];
    % Create line definition (domain 1).
    line_1 = [idx_domain_1, [idx_domain_1(2:end); idx_domain_1(1)]];
    n_line_1 = size(line_1, 1);
    
    % Get point indices in correct order (domain 2).
    pt_top_mid = find((point_list(:, 1) >= x0 & point_list(:, 1) <= (x0 + h)) &...
                      point_list(:, 2) == 0);
    pt_bot_mid = flipud(find((point_list(:, 1) >= x0 & point_list(:, 1) <= (x0 + h)) &...
                              point_list(:, 2) == bnd(4))...
                       );
    idx_domain_2 = [pt_top_mid; pt_bot_mid];
    % Create line definition (domain 2).
    line_2 = [idx_domain_2, [idx_domain_2(2:end); idx_domain_2(1)]];
    n_line_2 = size(line_2, 1);
    
    % Get point indices in correct order (domain 3).
    pt_top_right = find(point_list(:, 1) >= (x0 + h) & point_list(:, 2) == 0);
    pt_bot_right = flipud(find(point_list(:, 1) >= (x0 + h) & ...
                               point_list(:, 2) == bnd(4))...
                         );
    idx_domain_3 = [pt_top_right; pt_bot_right];
    % Create line definition (domain 3).
    line_3 = [idx_domain_3, [idx_domain_3(2:end); idx_domain_3(1)]];
    n_line_3 = size(line_3, 1);
    
    % Create line input ([n x 3] matrix) line-structure.
    line_list = [line_1; line_2; line_3];
    n_line = n_line_1 + n_line_2 + n_line_3;
    line_input = [(1:n_line).', line_list];
    
    %% Extract boundary line indices.
    
    pt_top = find(point_list(:, 2) == 0);
    pt_bot = find(point_list(:, 2) == bnd(4));
    pt_left = find(point_list(:, 1) == bnd(1));
    pt_right = find(point_list(:, 1) == bnd(2));
    
    idx_line_top = find(all(ismember(line_list, pt_top).'));
    idx_line_bot = find(all(ismember(line_list, pt_bot).'));
    idx_line_left = find(all(ismember(line_list, pt_left).'));
    idx_line_right = find(all(ismember(line_list, pt_right).'));
    
    %% Prepare physical entity handling.
    
    point_idx = ((n_point_bnd + n_point_dike + 1):n_point).';
    point_idx = find(ismember(sort_idx, point_idx)); % due to sorting
    tmp = cumsum([n_line_1, n_line_2, n_line_3]);
    line_loop_idx = {1:tmp(1), tmp(1)+1:tmp(2), tmp(2)+1:tmp(3)}.';
    
    %% Summarize output.
    
    info = struct();
    info.point_list = point_input;
    info.line_list = line_input;
    info.idx_ele_point = point_idx;
    info.idx_line_loop = line_loop_idx;
    info.dom_name = {'left', 'middle', 'right'};
    info.idx_bnd = {idx_line_top, idx_line_bot, idx_line_left, idx_line_right}.';
    info.bnd_name = {'top', 'bot', 'left', 'right'};
end

function [] = writeGEO(name, info)
    % Create Gmsh geometry (text) file and define higher order entities.
    
    % Fetch info.
    point = info.point_list;
    line = info.line_list;
    point_idx = info.idx_ele_point;
    loop_idx = info.idx_line_loop;
    bnd_idx = info.idx_bnd;
    dom_name = info.dom_name;
    bnd_name = info.bnd_name;
    
    n_point = size(point, 1);
    n_line = size(line, 1);
    
    fileID = fopen(name, 'w');
 
    %% Handle geometric entities.
    
    % Note: multiple defined lines are merged automatically in Gmsh!
    
    % Add points.
    for ii = 1:n_point
        fprintf(fileID, 'Point(%d) = {%d, %d, %d, %d};\n', ...
            point(ii,:));
    end
    
    % Add lines.
    fprintf(fileID, '\n');
    for ii = 1:n_line
        fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
            line(ii,:));
    end
        
    % Add line loops (polygonal chains).
    idx_start = n_line + 9;
    id_line_loop = zeros(length(loop_idx), 1);
    fprintf(fileID, '\n');
    for ii = 1:length(loop_idx)
        cur_n_line = size(loop_idx{ii}, 2);
        id_line_loop(ii) = idx_start + ii;
        fprintf(fileID, ...
            ['Line Loop(%d) = {', repmat('%d, ', 1, cur_n_line-1), '%d};\n'], ...
            id_line_loop(ii), loop_idx{ii});
    end
    
    % Add geomentric surface.
    idx_start = id_line_loop(end) + 9;
    id_plane_surface = zeros(length(loop_idx), 1);
    fprintf(fileID, '\n');
    for ii = 1:length(loop_idx)
        id_plane_surface(ii) = idx_start + ii;
        fprintf(fileID, ...
            'Plane Surface(%d) = {%d};\n', ...
            id_plane_surface(ii), id_line_loop(ii));
    end
    
    %% Handle physical entities.
    
    % Add physical points.
    fprintf(fileID, '\n');
    fprintf(fileID, ...
            ['Physical Point("ele") = {', ...
                repmat('%d, ', 1, length(point_idx)-1), '%d};\n'], ...
            point_idx);
    
    % Add physical lines.
    fprintf(fileID, '\n');
    for kk = 1:length(bnd_name)
        fprintf(fileID, ...
            ['Physical Line("%s") = {', ...
                repmat('%d, ', 1, length(bnd_idx{kk})-1), '%d};\n'], ...
            bnd_name{kk}, bnd_idx{kk});
    end
    
    % Add physical surfaces.
    fprintf(fileID, '\n');
    for ii = 1:length(loop_idx)
        fprintf(fileID, ...
            'Physical Surface("%s") = {%d};\n', ...
            dom_name{ii}, id_plane_surface(ii));
    end
    
    % Remove multiples.
    fprintf(fileID, '\n');
    fprintf(fileID, 'Coherence;');
    
    fclose(fileID);
end