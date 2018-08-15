function mesh = loadGmsh(name, verbosity)
    % Read out mesh information from Gmsh .msh structure.
    %
    % SYNTAX
    %   mesh = loadGmsh(name[, verbosity])
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file name (with extention) to import 
    %            from.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing: 
    %            dimension, 
    %            vertex list, vertex2simplex list, 
    %            simplex2parameter domain list, parameter vector,
    %            vertex2boundary egde list,
    %            boundary edge2boundary parts, boundary part names
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed
    %
    % REMARKS
    %   To be able to handle arbitrary .msh inputs it some requirements on
    %   the gmsh files are made:
    %
    %   - each physical_line which has a physical line name is considered
    %     to be a domain boundary 
    %     (other physical lines are ignored/treated as internal boundaries)
    %   - the names can be arbitrary char strings
    %     (these will be used to link boundary conditions see Fe.getBndDOF)
    %
    %   - the physical_surface names consist of characters denoting the 
    %     corresponding domain/cell conductivivies.
    %     (char strings of letters will be refused - only 'Inf' is allowed)
    %
    %   Note that the obtained boundary edge list does not refer to the
    %   complete list of edges within the given mesh!
    %   To achieve this relations applying Mesh.appendEdgeInfo() is
    %   absolutely necessary.
    
    %% Check input.
    
    if nargin < 2
       verbosity = false;
    else
       assert(islogical(verbosity), ...
           'verbosity - logical specifying the output verbosity, expected.'); 
    end
    
    if verbosity
       fprintf('read Gmsh output ... '); 
    end
    
    %% Load .msh file.
    
    % Get whole content from file.
    file_content = fileread(name);
    
    % Remove obsolete pattern.
    file_content = regexprep(file_content, '"', '');
    
    % Store each line of the input file in a seperate cell.
    file_content = textscan(file_content, '%s', 'delimiter', '\n', ...
        'whitespace', '');
    file_content = file_content{1};
    
    % Check file.
    if isempty(file_content)
        error('Empty file.');

    elseif ~strcmp(file_content{1}, '$MeshFormat')
        error(['File to import from seems not to be an .msh ', ...
            'ascii formatted file.']);
    end
    
    % Go through content and search for leading keywords
    % (see gmsh file format documentation) to separate information blocks.
    idx_node_start = find(cellfun(@(x) strcmp(x, '$Nodes'), ...
        file_content)) + 2;
    idx_node_stop = find(cellfun(@(x) strcmp(x, '$EndNodes'), ...
        file_content)) - 1;
    idx_element_start = find(cellfun(@(x) strcmp(x, '$Elements'), ...
        file_content)) + 2;
    idx_element_stop = find(cellfun(@(x) strcmp(x, '$EndElements'), ...
        file_content)) - 1;
    idx_phys_start = find(cellfun(@(x) strcmp(x, '$PhysicalNames'), ...
        file_content)) + 1;
    idx_phys_stop = find(cellfun(@(x) strcmp(x, '$EndPhysicalNames'), ...
        file_content)) - 1;
    
    % Get vertex information.
    % Information structure:
    % $
    % Number of vertices
    % vtx. number, x-coord, y-coord, z-coord
    % $End
    vtx_content = file_content(idx_node_start:idx_node_stop);
    
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
    % ele. type == 1 -> edge
    %           == 2 -> triangle
    %
    %    1. tag == physical entity id
    %                -> phys. line / phys. surface
    %    2. tag == elementary geometrical entity id
    %                -> straight line / plane surface id
    % (As two types of information are included in the element list, the
    % above trick doesn't work here.)
    ele_content = file_content(idx_element_start:idx_element_stop);
    ele_content = cellfun(@(x) {sscanf(sprintf('%s ', x), '%f').'}, ...
                        ele_content);
    ele_types = cellfun(@(x) x(2), ele_content);
    assert(all(ele_types <= 2), ...
        'Gmsh file contains 3D data (or unsupported 2D element format).');
    
    % Exclude boundary edge information.
    edge_content_idx = ele_types == 1;
    edge_content = cell2mat(ele_content(edge_content_idx));
    assert(~isempty(find(edge_content_idx, 1)), ...
        ['Edge information is missing. Please add (straight) line and ', ...
        'appropriate physical line information to gmsh file.']);
    
    % Exclude cell/simplex information.
    cell_content = cell2mat(ele_content(~edge_content_idx));
    parameter_domains = cell_content(:,4);
    
    % Exclude physical properties of cells and edges.
    if isempty(idx_phys_start)
        warning(['No physical properties for cells (edges) are set. ', ...
            'Set mesh.params = 1 and assuming homogeneous Neumann ', ...
            'boundary conditions on every given edge.']);
        params = ones(size(cell_content, 1), 1);
    else
        % Information structure:
        % $
        % Number of phys. names
        % phys. dimension, phys. number, phys. name
        % $End
        % phys. dimension == 1 -> edge
        %                 == 2 -> triangle
        % phys. number == physical entity id (-> see above)
        phys_content = file_content(idx_phys_start:idx_phys_stop);
        phys_num = str2double(phys_content{1});
        phys_content = cellfun(@strsplit, phys_content(2:end), ...
                          'UniformOutput', false);
        phys_content = vertcat(phys_content{:});
        
        % Separate phys. dimensions.
        phys_ids = str2double(phys_content(:,1)); 
        phys_edge_id = phys_ids == 1;
        phys_dom_id = phys_ids == 2;
        
        % Check consistencies.
        if length(unique(edge_content(:,4))) ~= length(find(phys_edge_id))
           warning(['The number of physical lines differs from the ', ...
               'number of physical line names. Lines which are not ', ...
               'associated with a physical line name are considered ', ...
               'as internal boundary and, hence, are ignored from now on.']); 
        end
        assert(phys_num == length(find([phys_edge_id, phys_dom_id])), ...
           ['The number of physical domains differs from the number', ...
           'of detected edges and surface domains.']); 
        assert(length(find(phys_dom_id)) == length(unique(parameter_domains)), ...
            ['Number of physical domains and domains in mesh do not ', ...
            'match. Make sure that each (plane) surface is associated ', ...
            'with a respective physical surface whose (physical) ', ...
            'name is the layer conductivity (given as string).']);
        
        % Exclude edge/bnd information.
        phys_edge_content = phys_content(phys_edge_id,:);
        phys_edge_id = str2double(phys_edge_content(:,2));
        edge2bnd_id = edge_content(:,4) == phys_edge_id.';
        phys_edge_names = phys_edge_content(:,3).';
        
        % Exclude domain information and set up cell/parameter information.
        phys_dom_content = str2double(phys_content(phys_dom_id,:));
        assert(all(all(~isnan(phys_dom_content))), ...
            ['Converting domain (strings) info into double failed. ', ...
            'Make sure that each character string of physical names ', ...
            'only contains the domain conductivity values.']);
        phys_dom_id = phys_dom_content(:,2);
        param2dom_id = parameter_domains == phys_dom_id.';
        params = param2dom_id * phys_dom_content(:,3);
    end
    
    %% Summarize.
    
    mesh = struct();
    mesh.dim = max(ele_types);
    mesh.vertices = vertices;
    mesh.cell2vtx = cell_content(:, [6, 7, 8]);
    mesh.parameter_domains = parameter_domains;
    mesh.params = params;
    mesh.bnd_edge2vtx = edge_content(:, [6, 7]);
    mesh.bnd_edge_part = edge2bnd_id;
    mesh.bnd_edge_part_name = phys_edge_names;
    
    % Check consistencies.
    assert(~isempty(mesh.bnd_edge2vtx), ...
        'Domain boundaries could not be assigned.');
    assert(~isempty(mesh.parameter_domains), ...
        'Parameter domain(s) could not be assigned.');
    
    if verbosity
       fprintf('done.\n'); 
    end
end