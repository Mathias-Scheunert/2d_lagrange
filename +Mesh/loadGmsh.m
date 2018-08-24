function mesh = loadGmsh(name, args)
    % Read out mesh information from Gmsh .msh file (ascii, format Ver. 2).
    %
    % Additionally, an uniform mesh refinement can be applied before data
    % import.
    %
    % The function supports two operation modes:
    %   1) .msh file includes physical entities:
    %       - only the entities which are additionally equipped with a
    %       physical name are stored. Other entities will be ignored.
    %
    %   2) .msh file includes no physical entities:
    %       - all occuring geometric entities will be stored.
    %
    % SYNTAX
    %   mesh = loadGmsh(name[, ref, verbosity])
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file name (with extention) to import 
    %            from.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing:
    %       1) dim                                      - scalar
    %          vertices, points, edges, faces, volumes
    %            dimension, 
    %            vertex list, vertex2simplex list, 
    %            simplex2parameter domain list, parameter vector,
    %            vertex2boundary egde list,
    %            boundary edge2boundary parts, boundary part names
    %
    % OPTIONAL PARAMETER
    %   args.verbosity ... Logical, denoting if current status should be
    %                      printed
    %   args.ref       ... Scalar, denoting the number of uniform 
    %                      refinement steps.
    %
    % REMARKS
    %
    %   If no physical entities at all are set for geometrical entities in
    %   mesh, Gmsh seems to store every occuring geometric element in the 
    %   element list of the .msh file (e.g. all points are listed).
    %   -> there physical entity id is set to '0' by default
    %
    %   However, if physical entities are set for a subset of geometric
    %   entities (e.g. some points in mesh) Gmsh seems not to mention the
    %   remaining gometric entities anymore!
    %
    %   To be able to handle arbitrary .msh inputs in order to use their
    %   information for setting up FEM simulations, some requirements 
    %   on the Gmsh files / meshes are made:
    %
    %   - physical entities including physical entity names for all
    %     relevant geometric entities have to be set!
    %
    %   - each physical_line has a physical line name. These are considered
    %     to be domain boundaries (2D mesh).
    %     (if occur: physical lines without name are ignored)
    %   - the names can be arbitrary char strings
    %
    %   - each physical_surface has a physical surface name. These are
    %     consired to be parameter domains (2D mesh) or boundaries 
    %     (3D mesh).
    %   - the names can be arbitrary char strings
    %
    %   - each physical volume has a physical volume name. These are
    %     considere to be parameter domains (3D mesh).
    %   - the names can be arbitrary char strings
    %
    %   Note 2D mesh:
    %   The obtained (boundary) edge list does only comprise a subset of
    %   the entire list of edges in the mesh.
    %   To achieve this relations additional routines 
    %   (e.g. Mesh.appendEdgeInfo) are required.
    %
    %   Note 3D mesh:
    %   The obtained (boundary) face list does only comprise a subset of
    %   the entire list of faces in the mesh.
    %   To achieve this relations additional routines are required.
    %
    %   Depending on the applied gmsh meshing variant during the mesh
    %   generation (i.e gmsh -2 [-> 2D mesh] or gmsh -3 [-> 3D]) the
    %   element list always comprises the complete list of triangels [2D] 
    %   or thedrahera [3D].
    %   The node list always comprises the complete list of vertices,
    %   forming the mesh.
    
    %% Check input.
    
%{
    assert(isstruct(args) && all(isfield(args, {'ref', 'verbosity'})));
    assert(ischar(name), ...
        ['name - Character of Gmsh file name to load ', ...
        '(including file extention!).']);
    assert(isscalar(args.ref) && ~islogical(args.ref) && args.ref >= 0, ...
        ['ref - Scalar, denoting the number of uniform ref steps, ', ...
        'expected.']);
    assert(islogical(args.verbosity), ...
        ['verbosity - logical, denoting if status should be printed, ', ...
        'expected']);
%}
    % Not required as already done in Mesh.initMesh().
    
    %% Refine uniformly.
    
    if args.verbosity
        fprintf('Refine mesh ... '); 
    end
    
    name_tmp = [name(1:end-4), '_tmp', name(end-3:end)];
    copyfile(name, name_tmp);
    gmsh_path = dir('**/gmsh');
    for i = 1:args.ref
        system([gmsh_path.folder, '/gmsh -refine -v 0 ', ...
            name_tmp]);
    end
    if args.verbosity
       fprintf('done.\n');
    end
    
    %% Load .msh file.
    
    if args.verbosity
       fprintf('Load Gmsh output ... '); 
    end
    
    % Get whole content from file.
    file_content = fileread(name_tmp);
    delete(name_tmp);
    
    % Remove obsolete pattern.
    file_content = regexprep(file_content, '"', '');
    
    % Store each line of the input file in a seperate cell.
    file_content = textscan(file_content, '%s', 'delimiter', '\n', ...
        'whitespace', '');
    file_content = file_content{1};
    
    % Check file.
    assert(~isempty(file_content), ...
        'Empty file detected.');
    
    % Go through content and search for leading keywords which separates
    % information blocks (see gmsh file format documentation).
    name_list = {'$MeshFormat', '$EndMeshFormat', ...
                 '$PhysicalNames', '$EndPhysicalNames', ...
                 '$Nodes', '$EndNodes', ...
                 '$Elements', '$EndElements'};
    [~, name_in_file] = ismember(file_content, name_list);
    [idx_in_file, ~, idx_in_name] = find(name_in_file);
    assert(all(ismember([1:2, 5:8], idx_in_name)), ...
        ['Could not observe relevant name tags in the provided ', ...
        'file. Check if input is a Gmsh .msh file (given in ASCII ', ...
        'format version 2']);
    if ~any(ismember(3:4, idx_in_name))
        warning(sprintf(['Given file does not contain physical names.', ...
            '\nExport mode 1) will be used (see function doc).']));
    end
    idx_format_start = idx_in_file(idx_in_name == 1) + 1;
    idx_format_stop = idx_in_file(idx_in_name == 2) - 1;
    idx_phys_start = idx_in_file(idx_in_name == 3) + 1;
    idx_phys_stop = idx_in_file(idx_in_name == 4) - 1;
    idx_node_start = idx_in_file(idx_in_name == 5) + 2;
    idx_node_stop = idx_in_file(idx_in_name == 6) - 1;
    idx_element_start = idx_in_file(idx_in_name == 7) + 2;
    idx_element_stop = idx_in_file(idx_in_name == 8) - 1;

    % Check if file type is supported.
    format_content = file_content(idx_format_start:idx_format_stop);
    format_content = strsplit(format_content{:});
    if sscanf(format_content{1}, '%d') ~= 2 ...
        % As in principal all keywords were found, import may work.
        warning(['File format differs from version 2. Import from ', ...
                'different versions may lead to wrong results.']);
    end
       
    %% Get element information.

    % Information structure:
    % $
    % Number of ele.
    % ele. number, ele. type, number of ele. tags, [tag list], [vertex list]
    % $End
    %
    % Note: 
    %   [tag list] as many columns/entries as given by number of ele. tags
    %   [vertex list] differs in length for edges and cells
    %
    % ele. type == 1  -> edge
    %           == 2  -> triangle
    %           == 4  -> tetrahedron
    %           == 15 -> point
    %
    %    1. tag == physical entity id
    %                -> phys. line / phys. surface / phys. volume
    %    2. tag == elementary geometrical entity id
    %                -> straight line / plane surface id / volume id
    ele_content_num = file_content(idx_element_start - 1);
    ele_content_num = str2double(ele_content_num{:});
    ele_content_tmp = file_content(idx_element_start:idx_element_stop);
    n_ele_cells = length(ele_content_tmp);
    ele_content = cell(n_ele_cells, 1);
    ele_content_tmp = char(ele_content_tmp);
    for ii = 1:n_ele_cells
       ele_content{ii} = sscanf(ele_content_tmp(ii,:), '%f').';
    end
    assert(size(ele_content, 1) == ele_content_num, ...
        'Mismatch between obtained and expected number of elements.');

    % Get all occuring element types and check if mesh type is supported.
    ele_types = cellfun(@(x) x(2), ele_content);
    supported_ele_types = [1, 2, 4, 15];
    found_ele_types = unique(ele_types);
    assert(all(ismember(found_ele_types, supported_ele_types)), ...
        sprintf(['File contains unsupported element types.\n', ...
        'Currently supported:', ...
        '\n id \t type', ...
        '\n 1 \t edge', ...
        '\n 2 \t triangle', ...
        '\n 4 \t tetrahedron', ...
        '\n 15 \t point']));
    assert(any(ismember(found_ele_types, [2, 4])), ...
        ['File solely contains point and/or edge information. ', ...
        'Only 2D meshes (containing triangle-elements) or 3D meshes ', ...
        '(containing thetrahedra-elements) are supported.']);
    
    % Determine the mesh dimension.
    if any(found_ele_types == 4)
        dim = 3;
    else
        dim = 2;
    end

    % Extract point information.
    point_content_idx = ele_types == 15;
    points = cell2mat(ele_content(point_content_idx));
    
    % Extract edge information.
    edge_content_idx = ele_types == 1;
    edges = cell2mat(ele_content(edge_content_idx));

    % Extract surface information.
    face_content_idx = ele_types == 2;
    faces = cell2mat(ele_content(face_content_idx));

    % Extract volume information.
    volume_content_idx = ele_types == 4;
    volumes = cell2mat(ele_content(volume_content_idx));

    % Summarzie.
    % Don't change order -> see paragraph "Get physical entities ..."
    element_name = {'point2vtx', 'edge2vtx', 'face2vtx', 'volume2vtx'};
    element_content = {points; edges; faces; volumes};
    req_cols = 1:4; % only these colums contain element2vtx relations
    element_empty = cellfun(@isempty, element_content);
    element_exist = find(~element_empty);
    element_n_tag = cellfun(@(x) {x(:,3)}, ...
                        element_content(~element_empty));
    assert(all(cellfun(@(x) all(x(1) == x), element_n_tag)), ...
        ['Number of tags within the list of one or more elements ', ...
        'types differs.']);

    %% Get vertex information.

    % Information structure:
    % $
    % Number of vertices
    % vtx. number, x-coord, y-coord, z-coord
    % $End
    vtx_content_num = file_content(idx_node_start - 1);
    vtx_content_num = str2double(vtx_content_num{:});
    vtx_content = file_content(idx_node_start:idx_node_stop);
    
    % Use a trick to speed up transforming the strings within the cells
    % into an array of numbers (as number of rows and cols don't change).
    cols = size(str2double(strsplit(vtx_content{1})), 2);
    rows = size(vtx_content, 1);
    vtx_content = sprintf('%s ', vtx_content{:});
    vertices = reshape(sscanf(vtx_content, '%f'), cols, rows).';
    vertices(:,1) = [];
    n_vtx = size(vertices, 1);
    assert(n_vtx == vtx_content_num, ...
        'Mismatch between obtained and expected number of vertices.');

    % Keep only relevant information.
    switch dim
        case 2
            % Obtain non-zero coordinate direction.
            eps_tol = pick(1, 5e0, 1e2);
            cur_dim = ~(sum(abs(vertices), 1) < eps * n_vtx * eps_tol);
            vertices = vertices(:, cur_dim);
        case 3
            % Nothing to be done.
    end

    %% Get physical entities of cells and edges and / or create output.

    if isempty(idx_phys_start)

        % Summarize.
        mesh = struct();
        mesh.dim = dim;
        mesh.vertices = vertices;
        for ii = 1:length(element_exist)
            mesh.(element_name{element_exist(ii)}) = ...
                element_content{element_exist(ii)}...
                    (:,end-(req_cols(element_exist(ii))-1):end);
        end
        switch dim
            case 2
                mesh.notes = ['Except of "vertices" and "faces", ', ...
                    'fields may not refer to the complete list of ', ...
                    'respective elements.'];
            case 3
                mesh.notes = ['Except of "vertices" and "volumes", ', ...
                    'fields may not refer to the complete list of ', ...
                    'respective elements.'];
        end
       
    else
        error('Not completely implemented yet.');

        % Information structure:
        % $
        % Number of phys. names
        % phys. dimension, phys. number, phys. name
        % $End
        % phys. dimension == 0 -> point
        %                 == 1 -> line (edge)
        %                 == 2 -> face (triangle)
        %                 == 3 -> volume (tetrahedron) 
        % phys. number == physical entity id (-> see above)
        phys_names = {'point', 'edge', 'face', 'volume'};
        phys_content = file_content(idx_phys_start:idx_phys_stop);
        phys_num = str2double(phys_content{1});
        phys_content = cellfun(@strsplit, phys_content(2:end), ...
                          'UniformOutput', false);
        phys_content = vertcat(phys_content{:});

        % Get all occuring physical dim types.
        % Don't change order
        phys_types = cellfun(@(x) str2double(x), phys_content(:,1));
        supported_phys_types = [0, 1, 2, 3];
        found_phys_types = unique(phys_types);
        phys_exist = ismember(supported_phys_types, found_phys_types);
        assert(all(ismember(found_phys_types, supported_phys_types)), ...
            sprintf(['File contains unsupported physical element types.\n', ...
            'Currently supported:', ...
            '\n id \t type', ...
            '\n 0 \t point', ...
            '\n 1 \t edge', ...
            '\n 2 \t triangle', ...
            '\n 3 \t tetrahedron']));

        % Separate phys. dimensions (=geometrical entity types).
        phys_point_id = phys_types == 0;
        phys_points = phys_content(phys_point_id,:);
        phys_edge_id = phys_types == 1;
        phys_edges = phys_content(phys_edge_id,:);
        phys_face_id = phys_types == 2;
        phys_faces = phys_content(phys_face_id,:);
        phys_volume_id = phys_types == 3;
        phys_volumes = phys_content(phys_volume_id,:);
        phys_ids = {phys_points; phys_edges; phys_faces; phys_volumes};
        assert(length(phys_types) == phys_num, ...
            ['Mismatch between obtained and expected number of ', ...
            'physical entities.']);

        % Check consistencies.
        phys_id_in_ele_cont = cell(4, 1);
        phys_id_in_ele_cont(element_exist) = ...
            cellfun(@(x) {unique(x(:,4))}, element_content(element_exist));
        pt_name_vs_pt_ele = cellfun(@(x, y) {unique(x(:,4)) == ...
                        unique(str2double(y(:,2)))}, ...
                            element_content(element_exist), ...
                            phys_ids(phys_exist));
        pt_mismatch = cellfun(@(x) {~all(find(x))}, pt_name_vs_pt_ele);
        if ~isequal(~element_empty(:), phys_exist(:))
            waring(['Mismatch between physical names and physical ', ...
                'entities observed.']);
        elseif ~all(vertcat(pt_name_vs_pt_ele{:}))
            waring(['Mismatch between physical entity ids within ', ...
                'physical names list and physical entity ids related ', ...
                'to geometric entity list.']);
        end
        
        % TODO: continue implementing.
        % -> check consistency of separate entity types.
        % -> only exclude those information which are related to physical
        % names

        if length(unique(edges(:,4))) ~= length(find(phys_edge_id))
           warning(sprintf(['\nThe number of physical lines differs ', ...
               'from the number of physical line names. \nLines which ', ...
               'are not associated with a physical line name are ', ...
               'considered as internal boundary and, hence, are ', ...
               'ignored from now on.']));

           % Reduce edge_content.
           obsolete_idx = ~ismember(edges(:,4), ...
               str2double(phys_content(phys_edge_id,2)));
           edges(obsolete_idx,:) = [];
        end
        assert(length(unique(faces(:,4))) == length(find(phys_volume_id)), ...
            ['Number of physical surfaces differs from the number of ', ...
            'physical surface names.']);

        % Exclude edge/bnd information.
        % (Provide a new progressive indexing, starting from 1 - that
        % matches the ordering of names.)
        phys_edge_content = phys_content(phys_edge_id,:);
        phys_edge_id = str2double(phys_edge_content(:,2));
        edge2bnd_id = edges(:,4) == phys_edge_id.';
        edge2bnd_id = edge2bnd_id * (1:length(phys_edge_id)).';
        phys_edge_names = phys_edge_content(:,3);

        % Exclude domain information and set up cell/parameter information.
        % (Provide a new progressive indexing.)
        phys_dom_content = phys_content(phys_volume_id,:);
        phys_volume_id = str2double(phys_dom_content(:,2));
        param2dom_id = faces(:,4) == phys_volume_id.';
        param2dom_id = param2dom_id * (1:length(phys_volume_id)).';
        parameter_domain_names = phys_dom_content(:,3);

        % Summarize.
        mesh = struct();
        mesh.dim = dim;
        mesh.vertices = vertices;
        mesh.cell2vtx = faces(:, [6, 7, 8]);
        mesh.parameter_domain = param2dom_id;
        mesh.parameter_domain_name = parameter_domain_names;
        mesh.bnd_edge2vtx = edges(:, [6, 7]);
        mesh.bnd_edge_part = edge2bnd_id;
        mesh.bnd_edge_part_name = phys_edge_names;
    end
    
    if args.verbosity
       fprintf('done.\n'); 
    end
end