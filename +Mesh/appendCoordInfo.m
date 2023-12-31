function mesh = appendCoordInfo(mesh)
    % Append mappings between mesh infos and its coordinates.
    %
    % SYNTAX
    %   mesh = appendCoordInfo(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            -> no coordinate information contained.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, mesh information appended by the mappings between
    %            cells to vtx coords and edges to vtx coords.
    %
    % REMARKS
    %
    %   cell2cord: n * [3 x 2]
    %   edge2cord: k * [2 x 2]
    %
    %   n = number of cells
    %   k = number of edges
    %   rows   = index of vertex
    %   colums = first x, second y coord

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

    % Mapping between edges and its coordinates.
    edge_list = mat2cell(mesh.edge2vtx.', 2, ones(n_edge, 1)).';
    edge2cord = cellfun(@(x) [mesh.vertices(x, 1),  mesh.vertices(x, 2)], ...
        edge_list, 'UniformOutput', false);

    %% Summarize infos.

    mesh.cell2cord = cell2cord;
    mesh.edge2cord = edge2cord;
end
