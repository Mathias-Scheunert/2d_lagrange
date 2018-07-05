function mesh = appendMappings(mesh)
    % Append mappings between mesh infos and its coordinates.
    % 
    % SYNTAX
    %   mesh = appendMappings(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing mesh information, i.e. coordinates of
    %            vertices and its relation to the triangles and edges.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, mesh information appended by the mappings between
    %            the vertices coordinates and the edges and triangles.

    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge information, expected.');
       
    %% Add mappings.
    
    % Get sizes.
    n_cell = size(mesh.cell2vtx, 1);
    n_edge = size(mesh.edge2vtx, 1);
    
    % Mapping between cells and its coordinates (coords of the vertices).
    cell_list = mat2cell(mesh.cell2vtx.', 3, ones(n_cell, 1)).';
    cell2cord = cellfun(@(x) [mesh.vertices(x, 1),  mesh.vertices(x, 2)], ...
        cell_list, 'UniformOutput', false);
    
    % Check consistency (for refinement).
    if isfield(mesh, 'cell2cord')
        assert(isequal(mesh.cell2cord, cell2cord), ...
            'Refinement produced inconsistent grid.');
    end
    
    % Mapping between edges and its coordinates.
    % TODO: proof if this might be in conflict with external mesh info
    % (e.g. field mesh.edge2vtx missing)
    edge_list = mat2cell(mesh.edge2vtx.', 2, ones(n_edge, 1)).';
    edge2cord = cellfun(@(x) [mesh.vertices(x, 1),  mesh.vertices(x, 2)], ...
        edge_list, 'UniformOutput', false);
    
    %% Summarize infos.
    
    mesh.cell2cord = cell2cord;
    mesh.edge2cord = edge2cord;
end
