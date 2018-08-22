function mesh = loadGmsh(name, args)
    % Read out mesh information from Gmsh .msh structure.
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
    
    % Not required as already done in Mesh.initMesh().
    
    %% Refine uniformly.
    
    if args.verbosity
        fprintf('Refine mesh ... '); 
    end
    
    name_tmp = [name(1:end-4), '_tmp', name(end-3:end)];
    copyfile(name, name_tmp);
    for i = 1:args.ref
        system([pwd, '/+Mesh/External/gmsh -refine -v 0 ', ...
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
    if isempty(file_content)
        error('Empty file.');

    elseif ~strcmp(file_content{1}, '$MeshFormat')
        error(['File to import from seems not to be an .msh ', ...
            'ascii formatted file.']);
    end
    
    % Go through content and search for leading keywords
    % (see gmsh file format documentation) to separate information blocks.
    name_list = {'$PhysicalNames', '$EndPhysicalNames', ...
                 '$Nodes', '$EndNodes', ...
                 '$Elements', '$EndElements'};
    [~, name_in_file] = ismember(file_content, name_list);
    [idx_in_file, ~, idx_in_name] = find(name_in_file);
    assert(all(ismember(1:6, idx_in_name)), ...
        ['Not all relevant name tags could be found in the provided ', ...
        'mesh file.']);
    idx_phys_start = idx_in_file(idx_in_name == 1) + 1;
    idx_phys_stop = idx_in_file(idx_in_name == 2) - 1;
    idx_node_start = idx_in_file(idx_in_name == 3) + 2;
    idx_node_stop = idx_in_file(idx_in_name == 4) - 1;
    idx_element_start = idx_in_file(idx_in_name == 5) + 2;
    idx_element_stop = idx_in_file(idx_in_name == 6) - 1;
    
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
    % above trick doesn't work here - the number of columns are different.)
    ele_content_tmp = file_content(idx_element_start:idx_element_stop);
    n_ele_cells = length(ele_content_tmp);
    ele_content = cell(n_ele_cells, 1);
    ele_content_tmp = char(ele_content_tmp);
    for ii = 1:n_ele_cells
       ele_content{ii} = sscanf(ele_content_tmp(ii,:), '%f').';
    end
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
        warning(sprintf(['\nNo physical properties for cells (edges) ', ...
            'are set. \nSet mesh.params = 1 and assuming homogeneous ', ...
            'Neumann boundary conditions on every given edge.']));
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
           warning(sprintf(['\nThe number of physical lines differs ', ...
               'from the number of physical line names. \nLines which ', ...
               'are not associated with a physical line name are ', ...
               'considered as internal boundary and, hence, are ', ...
               'ignored from now on.']));
           
           % Reduce edge_content.
           obsolete_idx = ~ismember(edge_content(:,4), ...
               str2double(phys_content(phys_edge_id,2)));
           edge_content(obsolete_idx,:) = [];
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
    
    if args.verbosity
       fprintf('done.\n'); 
    end
end