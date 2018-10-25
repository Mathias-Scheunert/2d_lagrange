function mesh = loadGmsh(name, varargin)
    % Read out mesh information from Gmsh .msh file (ascii, format Ver. 2).
    %
    % Additionally, an uniform mesh refinement can be applied before data
    % import.
    %
    % Supported Gmsh version: 3.x
    %
    % The function supports two operation modes:
    %   1) .msh file includes physical entities:
    %       - only the entities which are additionally equipped with a
    %       physical name are stored. Other entities will be ignored.
    %
    %   2) .msh file includes no physical entities:
    %       - all available geometric entities will be stored.
    %
    % SYNTAX
    %   mesh = loadGmsh(name[, varargin])
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file name (with extention) to import 
    %            from.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed [default = false].
    %   ref       ... Scalar, denoting the number of uniform 
    %                 refinement steps [default = 0].
    %   force     ... Logical, forcing operation mode 2) even if phys.
    %                 names are available [default = false].
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing:
    %       1) dimension,
    %          vertices, points, edges, faces, volumes
    %       2) dimension,
    %          vertices,
    %          points, point names,
    %          simplex2vertex list, 
    %          simplex2parameter domain list, parameter domain names,
    %          boundary2vertex list,
    %          boundary2boundary(id) parts, boundary part names
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
    
    assert(ischar(name), ...
        ['name - Character of Gmsh file name to load ', ...
        '(including file extention!).']);    % Define possible input keys and its properties checks.
    input_keys = {'ref', 'verbosity', 'force'};
    assertRef = @(x) assert(isscalar(x) && ~islogical(x) && x >= 0, ...
        'ref - Scalar, denoting the number of uniform ref steps, expected.');
    assertLogic = @(x) assert(islogical(x), ...
        'verbosity - logical, denoting if status should be printed, expected');
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, 0, assertRef);
    parser_obj.addParameter(input_keys{2}, false, assertLogic);
    parser_obj.addParameter(input_keys{3}, false, assertLogic);
   
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;
    
    %% Check Gmsh version.
    
    gmsh_path = dir('**/gmsh');
    [~, version] = system([gmsh_path.folder, '/gmsh -version']);
    version = textscan(version, '%s', 'delimiter', '\n', ...
                       'whitespace', '');
    version = version{1}{end};
    assert(strcmp(version(1), '3'), ...
        'Wrong Gmsh version detected - please use version 3.x.');
    
    %% Refine uniformly.
    
    if args.verbosity
        fprintf('Refine mesh ... '); 
    end
    name_tmp = [name(1:end-4), '_tmp', name(end-3:end)];
    copyfile(name, name_tmp);
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
    % -> Leave elements of element_content untouched, only modify summary.
    element_content = {points; edges; faces; volumes};
    % -> Don't change order, see paragraph "Get physical entities ..."
    element_name = {'point2vtx', 'edge2vtx', 'face2vtx', 'volume2vtx'};
    req_cols = 1:4; % only these colums contain element2vtx relations
    element_exist = ~cellfun(@isempty, element_content);
    element_type_exist = find(element_exist);
    element_n_tag = cellfun(@(x) {x(:,3)}, ...
                        element_content(element_exist));
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

    if isempty(idx_phys_start) || args.force

        % Summarize.
        mesh = struct();
        mesh.dim = dim;
        mesh.vertices = vertices;
        for ii = 1:length(element_type_exist)
            mesh.(element_name{element_type_exist(ii)}) = ...
                element_content{element_type_exist(ii)}...
                    (:,end-(req_cols(element_type_exist(ii))-1):end);
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
        % Information structure of file header:
        % $
        % Number of phys. names
        % phys. dimension, phys. number, phys. name
        % $End
        % phys. dimension == 0 -> point
        %                 == 1 -> line (edge)
        %                 == 2 -> face (triangle)
        %                 == 3 -> volume (tetrahedron) 
        % phys. number == physical entity id (-> see above)
        % Note: Only those four dimension, i.e. entity types are considered!
        % -> Don't change order of the following definitions:
        supported_phys_types = [0, 1, 2, 3];
        supported_phys_type_names = {'point', 'edge', 'face', 'volume'};
        supported_phys_types_num = length(supported_phys_types);
        
        % Get physical entities content from file header.
        phys_content = file_content(idx_phys_start:idx_phys_stop);
        phys_num = str2double(phys_content{1});
        phys_content = cellfun(@strsplit, phys_content(2:end), ...
                          'UniformOutput', false);
        phys_content = vertcat(phys_content{:});

        % Get all occuring physical entity types.
        phys_types = cellfun(@(x) str2double(x), phys_content(:,1));
        assert(length(phys_types) == phys_num, ...
            ['Mismatch between obtained and expected number of ', ...
            'physical entities within header.']);
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

        % Separate phys. dimensions/entities (=geometrical entity types).
        phys_point_id = phys_types == 0;
        phys_points = phys_content(phys_point_id,:);
        phys_edge_id = phys_types == 1;
        phys_edges = phys_content(phys_edge_id,:);
        phys_face_id = phys_types == 2;
        phys_faces = phys_content(phys_face_id,:);
        phys_volume_id = phys_types == 3;
        phys_volumes = phys_content(phys_volume_id,:);
        phys_id_info = {phys_points; phys_edges; phys_faces; phys_volumes};
        phys_ids = cell(supported_phys_types_num, 1);
        phys_ids(phys_exist) = cellfun(@(x) {unique(str2double(x(:,2)))}, ...
                                     phys_id_info(phys_exist));

        % Check consistencies between names in header and names in element
        % list.
        phys_ids_in_ele_cont = cell(supported_phys_types_num, 1);
        phys_ids_in_ele_cont(element_type_exist) = ...
            cellfun(@(x) {unique(x(:,4))}, element_content(element_type_exist));
        
        % Check if complete entity types are not related to phys. names.
        if ~isequal(element_exist(:), phys_exist(:))
            tmp_names = supported_phys_type_names(...
                            element_exist(:) ~= phys_exist(:)...
                        );
            warning(sprintf(['Not every supported physical entity type is ', ...
                'related to physical names. Ignoring unrelated ', ...
                'entity types: \n', ...
                repmat('- %s\n', 1, length(tmp_names))], ...
                tmp_names{:}));
            
            % Remove obsolete entity types.
            element_exist(:) = element_exist(:) & phys_exist(:);
            element_type_exist = find(element_exist);
            element_content(~element_exist) = {[]};
            phys_ids_in_ele_cont(~element_exist) = {[]};
        end
        
        % Check if parts of existing entities are not related to
        % appropriate phys. names.
        pid_vs_pid_ele = cellfun(@(x, y) {ismember(x, y)}, ...
            phys_ids_in_ele_cont, phys_ids);
        pid_mismatch = cellfun(@(x) {~all(x)}, pid_vs_pid_ele);
        
        % Reduce information from element list to only comprise entity data
        % related to the observed physical names.
        if pid_mismatch{1}
            warning(['Number of physical point names differ from the ', ...
            'number of physical point ids. Ignoring unrelated point ids.']);
        
           % Reduce point_content.
           obsolete_idx = ~ismember(points(:,4), ...
               str2double(phys_content(phys_point_id, 2)));
           element_content{1}(obsolete_idx,:) = [];
        end
        if pid_mismatch{2}
           warning(['The number of physical line names differ from ', ...
               'the number of physical line ids. Ignoring unrelated ', ...
               'line ids.']);

           % Reduce edge_content.
           obsolete_idx = ~ismember(edges(:,4), ...
               str2double(phys_content(phys_edge_id, 2)));
           element_content{2}(obsolete_idx,:) = [];
        end
        if pid_mismatch{3}
            warning(['The number of physical face names differ from ', ...
               'the number of physical face ids. Ignoring unrelated ', ...
               'face ids.']);
        
            % Reduce face_content.
           obsolete_idx = ~ismember(faces(:,4), ...
               str2double(phys_content(phys_face_id, 2)));
           element_content{3}(obsolete_idx,:) = [];
        end
        if pid_mismatch{4}
            warning(['The number of physical volume names differ from ', ...
               'the number of physical volume ids. Ignoring unrelated ', ...
               'volume ids.']);
        
            % Reduce volume_content.
           obsolete_idx = ~ismember(volumes(:,4), ...
               str2double(phys_content(phys_volume_id, 2)));
           element_content{4}(obsolete_idx,:) = [];
        end

        % Generate entity/entity name mappings.
        % (Provide a new progressive indexing, starting from 1 - that
        % matches the ordering of names.)
        
        % Map point information.
        if ~isempty(element_content{1})
            phys_point_content = phys_content(phys_point_id,:);
            phys_point_id = phys_ids{1}(:);
            point_id_map = element_content{1}(:, 4) == phys_point_id.';
            point_id_map = point_id_map * (1:length(phys_point_id)).';
            phys_point_names = phys_point_content(:, 3);
        else
            point_id_map = [];
            phys_point_names = {[]};
        end
        
        % Map edge(= bnd in 2D) information.
        if ~isempty(element_content{2})
            phys_edge_content = phys_content(phys_edge_id,:);
            phys_edge_id = phys_ids{2}(:);
            edge_id_map = element_content{2}(:, 4) == phys_edge_id.';
            edge_id_map = edge_id_map * (1:length(phys_edge_id)).';
            phys_edge_names = phys_edge_content(:, 3);
        else
            edge_id_map = [];
            phys_edge_names = {[]};
            if dim == 2
                warning('No boundary information provided.');
            end
        end

        % Map face(= bnd in 3D or cell/parameter in 2D) information.
        if ~isempty(element_content{3})
            phys_face_content = phys_content(phys_face_id,:);
            phys_face_id = phys_ids{3}(:);
            face_id_map = element_content{3}(:, 4) == phys_face_id.';
            face_id_map = face_id_map * (1:length(phys_face_id)).';
            phys_face_names = phys_face_content(:, 3);
        else
            face_id_map = [];
            phys_face_names = {[]};
            if dim == 2
                warning('No domain information provided.');
            elseif dim == 3
                warning('No boundary information provided.');
            end
        end
        
        % Map volume information(= cell/parameter in 3D).
        if ~isempty(element_content{4})
            phys_volume_content = phys_content(phys_volume_id,:);
            phys_volume_id = phys_ids{4}(:);
            volume_id_map = element_content{4}(:, 4) == phys_volume_id.';
            volume_id_map = volume_id_map * (1:length(phys_volume_id)).';
            phys_volume_names = phys_volume_content(:, 3);
        else
            volume_id_map = [];
            phys_volume_names = {[]};
            if dim == 3
                warning('No domain information provided.');
            end
        end

        % Summarize.
        % Note: The name definitions are adapted to the 2d-lagrange project.
        mesh = struct();
        mesh.dim = dim;
        mesh.vertices = vertices;
        if ~isempty(point_id_map)
            mesh.point = point_id_map;
            mesh.point_names = phys_point_names;
        end
        switch dim
            case 2
                element_name = {'point2vtx', 'bnd_edge2vtx', ...
                                'cell2vtx', 'volume2vtx'};
                mesh.parameter_domain = face_id_map;
                mesh.parameter_domain_name = phys_face_names;
                mesh.bnd_edge_part = edge_id_map;
                mesh.bnd_edge_part_name = phys_edge_names;
            case 3
                element_name = {'point2vtx', 'edge2vtx', ...
                                'bnd_face2vtx', 'cell2vtx'};
                mesh.parameter_domain = volume_id_map;
                mesh.parameter_domain_name = phys_volume_names;
                mesh.bnd_face_part = face_id_map;
                mesh.bnd_face_part_name = phys_face_names;
        end
        for ii = 1:length(element_type_exist)
            mesh.(element_name{element_type_exist(ii)}) = ...
                element_content{element_type_exist(ii)}...
                    (:,end-(req_cols(element_type_exist(ii))-1):end);
        end
    end
   
    if args.verbosity
       fprintf('done.\n'); 
    end
end