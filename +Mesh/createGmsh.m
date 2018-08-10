function mesh = createGmsh(bnd, args)
    % Stores mesh information of a mesh created by Gmsh.
    %
    % Therefore, 
    %   1) a .geo input file is created 
    %   2) a .msh output file is obtained from running Gmsh
    %   3) information from .msh is converted
    %   4) .geo and .msh files are deleted
    %
    % INPUT PARAMETER
    %   bnd  ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   args ... Struct, may containing further vectorized (TX, Rx, topo)
    %            and scalar (ref) information.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing vertex2coordinates (vertices), 
    %            simplex2vertex (cell2vtx), boundary egde2coordinates 
    %            (bnd_edge2vtx), boundary egde 2 model domain boundary
    %            (bnd_edge_...) and simplex parameter domains 
    %            (parameter_domains).
    
    %% Check input.
    
    % Not required, as already done in Mesh.initMesh()
    
    %% Create Gmsh input file.
    
    gmsh_file = 'tmp_mesh';
    createGmshInput([gmsh_file, '.geo'], bnd, ...
        args.TX, args.RX, args.topo, args.verbosity);
    
    %% Run Gmsh.
    
    system([pwd, '/+Mesh/External/gmsh -2 ', gmsh_file, '.geo']);
    if args.verbosity
        fprintf('Refine mesh ... '); 
    end
    for i = 1:args.ref
        system([pwd, '/+Mesh/External/gmsh -refine -v 0 ', ...
            gmsh_file, '.msh']);
    end
    if args.verbosity
       fprintf('done.\n');
    end

    %% Import mesh information from .msh file.

    mesh = getGmshOutput([gmsh_file, '.msh'], ...
        args.verbosity);

    %% Clean up.
    
    delete([gmsh_file, '.geo'], [gmsh_file, '.msh']);
end

function createGmshInput(name, bnd, TX, RX, topo, verbosity)
    % Creates an .geo input file for Gmesh
    % 
    % SYNTAX
    %   createGmshInput(name, bnd, TX, RX, topo)
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file name to be created.
    %   bnd  ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   TX   ... Vector [n x 2], denoting the source position(s).
    %   RX   ... Vector [m x 2], denoting the receiver position(s).
    %   topo ... Vector [k x 2], denoting the descrete topography.
    %   ref  ... Scalar, denoting the number of uniform refinenments.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed
    %
    % REMARKS
    %  Note that the 2D grid will be referred to gmesh x and z coordinate:
    %  0,0 ------> 
    %      |      x = gmesh-x
    %      |
    %      |
    %      v
    %       y = gmesh-z
    %
    %  Note that gmsh can't handle point which lie on the boundary!
    %  I.g. if a TX/RX point is embedded in the surface but lies somewhere
    %  on a bounding line of that surface, gmsh will not be able to find a
    %  mesh for that situation.
    %  This might be due to the behavior, that gmsh meshes starting from
    %  lines to surfaces to volumes an can't handle additional points in
    %  surface which also occur on the already meshed lines.
    %
    %  To fix that issue, TX/RX points wich are not horizontaly aligned
    %  with a topography point, are expected to lie ON the boundary, such
    %  that the domain boundary will be adapted.
    
    if verbosity
       fprintf('\ndefine Gmesh input ... '); 
    end
        
    %% Set up domain basic geometry entities (points, lines).
        
    % Set domain boundaries
    X = [bnd(1:2), fliplr(bnd(1:2))]; 
    Y = [bnd(3:4); bnd(3:4)];
    Y = Y(:).';
    
    % Insert topography points.
    if ~isempty(topo)
        assert(all(topo(:,1) > bnd(1)) && all(topo(:,1) < bnd(2)), ...
            ['Points describing topography exceed the horizontal ', ...
            'domain extension.']);
    end
    
    % Set maximal element size.
    %                  geometry based              fix
    h_domain = pick(1, round(diff(bnd(1:2)) / 25), 10);
    
    % Create point input ([n x 5] matrix) structure for domain boundaries.
    n_domain = length(X);
    point_domain = zeros(n_domain, 5);
    for i = 1:n_domain
       point_domain(i, 1) = i ;       % point number
       point_domain(i, 2) = X(i);     % gmesh-x cooridnate.
       point_domain(i, 3) = 0;        % gmesh-y cooridnate.
       point_domain(i, 4) = Y(i);     % gmesh-z cooridnate.
       point_domain(i, 5) = h_domain; % element size around current point
    end
    
    % Add topography.
    if ~isempty(topo)
        % Sort points w.r.t. x-coordinate.
        topo = sortrows(topo);
        
        % Set maximal element size.
        % Use min point offset in x direction (~horizontal) & avoid zero.
        h_topo = max(min(diff(topo(:,1))), 0.1);
        
        % Create point input ([n x 5] matrix) structure for topography.
        n_topo = size(topo, 1);
        n_domain = n_topo + n_domain;
        point_topo = [(1:n_topo).', topo(:,1), zeros(n_topo, 1), ...
            topo(:,2), h_topo + zeros(n_topo, 1)];
        
        % Insert points into domain point list.
        point_domain = [point_domain(1, :); ...
                        point_topo; ...
                        point_domain(2:end, :)];
        point_domain(:,1) = 1:n_domain;
    end
    
    % Add TX/RX positions.
    if ~isempty([TX; RX])
        % Initialize.
        %                geometry based               fix
        h_TXRX = pick(1, round(diff(bnd(1:2)) / 500), 0.1);
        TXRX = [TX; RX];
        
        % Just reduce h, if points are already part of point_domain.
        [txrx_in_pd, idx_txrx_in_pd] = ismember(... 
            point_domain(:, [2, 4]), TXRX, 'rows');
        point_domain(txrx_in_pd,5) = h_TXRX;
        idx_txrx_not_in_pd = ~ismember(...
            1:size(TXRX, 1), idx_txrx_in_pd(idx_txrx_in_pd ~= 0));
        
        if any(idx_txrx_not_in_pd)
            % Check if left TX/RX points are horizontally aligned with
            % topography information. 
            if ~isempty(topo) && any(ismember(...
                    TXRX(idx_txrx_not_in_pd, 1), ...
                    topo(idx_txrx_not_in_pd, 1)))
               
                % Get point index.
                idx_internal = find(ismember(...
                                  TXRX(idx_txrx_not_in_pd, 1), ...
                                  topo(idx_txrx_not_in_pd, 1))...
                                );
                % Check, if point lies above the surface.
                if any(TXRX(idx_internal, 1) < topo(idx_internal, 1))
                    warning('TX/RX-Point lies outside the model domain.');
                end
                
                % Add points where finer mesh is needed to interior of 
                % the model domain.
                n_TXRX = length(idx_internal);
                point_internal = [n_domain + (1:n_TXRX).', ...
                    TXRX(idx_internal,1), ...
                    zeros(n_TXRX, 1), ...
                    TXRX(idx_internal,2), ...
                    h_TXRX + zeros(n_TXRX, 1)];
                
                % Remove point(s) from index vector.
                idx_txrx_not_in_pd(idx_internal) = false;
                
            else
                % Set points where finer mesh is needed to be part of the 
                % boundary of the model domain.
                
                % Insert points into domain point list.
                idx_txrx_on_domain = find(idx_txrx_not_in_pd);
                n_txrx_on_domain = length(idx_txrx_on_domain);
                n_domain = n_domain + n_txrx_on_domain;
                point_txrx_domain = [(1:n_txrx_on_domain).', ...
                                      TXRX(idx_txrx_not_in_pd, 1), ...
                                      zeros(n_txrx_on_domain, 1), ...
                                      TXRX(idx_txrx_not_in_pd, 2), ...
                                      h_TXRX + zeros(n_txrx_on_domain, 1)];
                point_tmp = point_domain(end,:);
                point_domain = [point_domain(1:end-1,:); point_txrx_domain];
                
                % Make sure, that points will form a closed point chain
                % (w.r.t. x-coo) which forms the domain boundary outline.
                [~, idx_sort_domain] = sort(point_domain(:,2));
                point_domain = point_domain(idx_sort_domain,:);
                point_domain = [point_domain; point_tmp];
                point_domain(:,1) = 1:n_domain;
                
                % Remove point(s) from index vector.
                idx_txrx_not_in_pd(idx_txrx_on_domain) = false;
            end            
        end
    end
    
    % Create line input ([n x 3] matrix) structure for domain boundaries.
    n_line = size(point_domain, 1);
    idx_start = 1:n_line;
    idx_end = [2:(n_line), 1];
    line = [1:n_line; idx_start; idx_end].';
    
    %% Create Gmsh geometry (text) file and define higher order entities.
    
    fileID = fopen(name, 'w');
    
    % Add geometry-describing point ids.
    for i = 1:n_domain
        fprintf(fileID, 'Point(%d) = {%d, %d, %d, %d};\n', ...
            point_domain(i,:));
    end
    
    % Add line ids.
    fprintf(fileID, '\n');
    for i = 1:n_line
        fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
            line(i,:));
    end
        
    % Define polygonal chain describing the complete domain boundary.
    id_line_loop = n_line + 10;
    fprintf(fileID, '\n');
    fprintf(fileID, ...
        ['Line Loop(%d) = {', repmat('%d, ', 1, n_line-1), '%d};\n'], ...
        id_line_loop, 1:n_line);
    
    % Define interior of this polygonal chain as closed volume (=domain).
    id_plane_surface = 10^(ceil(log10(id_line_loop)));
    fprintf(fileID, '\n');
    fprintf(fileID, ...
        'Plane Surface(%d) = {%d};\n', id_plane_surface, id_line_loop);
    
    % Add detached points.
    if any(idx_txrx_not_in_pd)
        fprintf(fileID, '\n');
        for i = 1:n_TXRX
            fprintf(fileID, 'Point(%d) = {%d, %d, %d, %d};\n', ...
                point_internal(i,:));
        end
    end
    
    % Add detached points to surface.
    % (Only in this case they will be considered by the meshing routine.)
    if any(idx_txrx_not_in_pd) && n_TXRX > 1
        fprintf(fileID, '\n');
        fprintf(fileID, ...
            ['Point{', repmat('%d, ', 1, n_TXRX-1),'%d} In Surface{%d};\n'], ...
            point_internal(:, 1), id_plane_surface);
    
    elseif any(idx_txrx_not_in_pd) && n_TXRX == 1
        fprintf(fileID, '\n');
        fprintf(fileID, ...
            'Point{%d} In Surface{%d};\n', ...
            point_internal(1, 1), id_plane_surface);
    end
    
    % Add boundary (physical) line ids.
%     id_phys_line = 2*10^(ceil(log10(id_plane_surface)));
    id_phys_line = 1;
    % xmin (left)
    fprintf(fileID, '\n');
    fprintf(fileID, 'Physical Line(%d) = {%d};\n', ...
        id_phys_line, line(end, 1));
    % xmax (right)
    fprintf(fileID, 'Physical Line(%d) = {%d};\n', ...
        id_phys_line+1, line(end-2, 1));
    % ymin (surface)
    if n_line > 4
        fprintf(fileID, ['Physical Line(%d) = {', ...
            repmat('%d, ', 1, n_line-4), '%d};\n'], ...
            id_phys_line+2, line(1:(end-3), 1));        
    elseif n_line == 4
        fprintf(fileID, 'Physical Line(%d) = {%d};\n', ...
            id_phys_line+2, line(1, 1));
    end
    % ymax (bottom)
    fprintf(fileID, 'Physical Line(%d) = {%d};\n', ...
        id_phys_line+3, line(end-1, 1));
    
    % Add physical surface id.
    id_phys_surface = 1;
    fprintf(fileID, '\n');
    fprintf(fileID, ...
        'Physical Surface(%d) = {%d};\n', id_phys_surface, id_plane_surface);
    
    fclose(fileID);
    
    if verbosity
       fprintf('done.\n'); 
    end
end

function mesh = getGmshOutput(name, verbosity)
    % Read out mesh information from Gmsh .msh structure.
    %
    % SYNTAX
    %   mesh = getGmshOutput(name[, verbosity])
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file (with extention) name to inport 
    %            from.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing vertex, vertex2simplex, boundary eges,
    %            simplex parameter domains.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed
    %
    % REMARKS
    %   To be able to handle arbitrary .msh inputs it is required
    %   that the physical_line tag for identifying domain boundaries are 
    %   set as follows:
    %       1 - boundary at xmin (left)
    %       2 - boundary at xmax (right)
    %       3 - boundary at ymin (top/surface)
    %       4 - boundary at ymax (bottom)
    %   -> other tags will be ignored / expected to be internal boundaries.
    %
    %   that the physical_surface tag for identifying areas of the equal 
    %   cell parameters are set as scalar(s).
    
    if verbosity
       fprintf('read Gmesh output ... '); 
    end
    
    %% Load .msh file.
    
    file_ID = fopen(name, 'r');
    
    % Store each line of the input file in a seperate cell.
    file_content = textscan(file_ID, '%s', 'delimiter', '\n', ...
        'whitespace', '');
    file_content = file_content{1};
    fclose(file_ID);
    
    % Check file.
    if isempty(file_content)
        error('Empty file.');

    elseif ~strcmp(file_content{1}, '$MeshFormat')
        error(['File to import from seems not to be an .msh ', ...
            'ascii formatted file.']);
    end
    
    % Go through content and search for leading keywords
    % (see gmsh file format documentation) to separate information blocks.
    idx_node_start = find(cellfun(@(x) strcmp(x,  '$Nodes'), ...
        file_content)) + 2;
    idx_node_stop = find(cellfun(@(x) strcmp(x,  '$EndNodes'), ...
        file_content)) - 1;
    idx_element_start = find(cellfun(@(x) strcmp(x,  '$Elements'), ...
        file_content)) + 2;
    idx_element_stop = find(cellfun(@(x) strcmp(x,  '$EndElements'), ...
        file_content)) - 1;
    
    % Get vertex information.
    vtx_content = file_content(idx_node_start:idx_node_stop);
%     vtx_content = cellfun(@(x) str2double(strsplit(x)), ...
%         vtx_content, 'UniformOutput', false);
%     vertices = cell2mat(vtx_content);
    
    % Use a trick to speed up transforming the strings within the cells
    % into an array of numbers,
    cols = size(str2double(strsplit(vtx_content{1})), 2);
    rows = size(vtx_content, 1);
    vtx_content = sprintf('%s ', vtx_content{:});
    vertices = reshape(sscanf(vtx_content, '%f'), cols, rows).';

    % Keep only relevant information.
    vertices = vertices(:, [2, 4]);
    
    % Consistency check.
    assert(size(vertices, 1) == ...
        str2double(cell2mat(file_content(idx_node_start - 1))), ...
        'Mismatch between obtained and expected number of vertices.');
       
    % Get element information.
    % (As two types of information are included in the element list, the
    % above trick doesn't work here.)
    ele_content = file_content(idx_element_start:idx_element_stop);
    ele_content = cellfun(@(x) {sscanf(sprintf('%s ', x), '%f').'}, ...
                        ele_content);
    
    % Exclude boundary edge information.
    bnd_edge_content_idx = cellfun(@(x) x(2) == 1, ele_content);
    bnd_edge_content = cell2mat(ele_content(bnd_edge_content_idx));
    
    % Exclude cell/simplex information.
    cell_content = cell2mat(ele_content(~bnd_edge_content_idx));
    
    %% Summarize.
    
    mesh = struct();
    mesh.dim = 2;
    mesh.vertices = vertices;
    mesh.cell2vtx = cell_content(:, [6, 7, 8]);
    mesh.parameter_domains = cell_content(:,4);
    mesh.bnd_edge2vtx = bnd_edge_content(:, [6, 7]);
    mesh.bnd_edge_xmin = bnd_edge_content(:,4) == 1;
    mesh.bnd_edge_xmax = bnd_edge_content(:,4) == 2;
    mesh.bnd_edge_ymin = bnd_edge_content(:,4) == 3;
    mesh.bnd_edge_ymax = bnd_edge_content(:,4) == 4;
    
    % Check consistencies.
    assert(~any([isempty(mesh.bnd_edge_xmin), ...
        isempty(mesh.bnd_edge_xmin), ...
        isempty(mesh.bnd_edge_xmin), ...
        isempty(mesh.bnd_edge_xmin)]), ...
        'Domain boundaries could not be assigned.');
    assert(~isempty(mesh.parameter_domains), ...
        'Parameter domain(s) could not be assigned.');
    
    if verbosity
       fprintf('done.\n'); 
    end
end