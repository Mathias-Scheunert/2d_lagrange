function mesh = loadGmsh(name, verbosity)
    % Read out mesh information from Gmsh .msh structure.
    %
    % SYNTAX
    %   mesh = loadGmsh(name[, verbosity])
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
    %   that the physical_surface names denote the corresponding cell
    %   conductivivies.
    
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
    idx_node_start = find(cellfun(@(x) strcmp(x,  '$Nodes'), ...
        file_content)) + 2;
    idx_node_stop = find(cellfun(@(x) strcmp(x,  '$EndNodes'), ...
        file_content)) - 1;
    idx_element_start = find(cellfun(@(x) strcmp(x,  '$Elements'), ...
        file_content)) + 2;
    idx_element_stop = find(cellfun(@(x) strcmp(x,  '$EndElements'), ...
        file_content)) - 1;
    idx_phys_dom_start = find(cellfun(@(x) strcmp(x,  '$PhysicalNames'), ...
        file_content)) + 1;
    idx_phys_dom_stop = find(cellfun(@(x) strcmp(x,  '$EndPhysicalNames'), ...
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
    assert(~isempty(find(bnd_edge_content_idx, 1)), ...
        'Boundary edge information is missing.');
    
    % Exclude cell/simplex information.
    cell_content = cell2mat(ele_content(~bnd_edge_content_idx));
    parameter_domains = cell_content(:,4);
    
    % Exclude physical properties of domains.
    if isempty(idx_phys_dom_start)
        warning(['No physical properties for domains (edges) are set. ', ...
            'Set params = 1.']);
        params = ones(size(cell_content, 1), 1);
    else
        % Get physical domain information.
        phys_dom_content = file_content(idx_phys_dom_start:idx_phys_dom_stop);
        phys_dom_num = str2double(phys_dom_content{1});
        phys_dom_content = cell2mat(cellfun(@(x) {sscanf(sprintf('%s ', x), '%f').'}, ...
                        phys_dom_content(2:end)));
                    
        % Check consistency.
        assert(phys_dom_num == length(unique(parameter_domains)), ...
            'Number of physical domains and domains in mesh do not match.');
        % TODO: check what happens of only surface and no resp. physical
        % surface are defined.
        param_merge = parameter_domains == phys_dom_content(:,2).';
        params = param_merge * phys_dom_content(:,3);
    end
    
    %% Summarize.
    
    mesh = struct();
    mesh.dim = 2;
    mesh.vertices = vertices;
    mesh.cell2vtx = cell_content(:, [6, 7, 8]);
    mesh.parameter_domains = parameter_domains;
    mesh.params = params;
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