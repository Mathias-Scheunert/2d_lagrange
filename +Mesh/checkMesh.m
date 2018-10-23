function [] = checkMesh(mesh)
    % Check that all elements of mesh are ordered (lexicographic).
    %
    % The lexicographic / ascending ordering of vertices is required for 
    % consistent transition of basis functions (especially for 
    % Raviart-Thomas elements) for two adjacent cells in mesh.
    %
    % SYNTAX
    %   [] = checkMesh(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    
    %% Check input.
    
    to_check = {'vertices', 'cell2vtx'};
    opt_to_check = {'edge2vtx', 'cell2edg'};
    assert(isstruct(mesh) && all(isfield(mesh, to_check)), ...
        'mesh - Struct, containing vertex and cell2vertex info, expected.');
    
    % Add further checks if possible.
    to_check = [to_check, opt_to_check(isfield(mesh, opt_to_check))];
    
    %% Check ordering of entities.
    
    
    fprintf(sprintf(['Check sorting of:\n', ...
        repmat('- %s\n', 1, length(to_check))], to_check{:}));
    
    vtx_sort = issortedrows(mesh.vertices, [2, 1]);
    cell2vtx_sort = issortedrows(mesh.cell2vtx.');
    if ismember('edge2vtx', to_check)
        edge2vtx_sort = issortedrows(mesh.edge2vtx.');
    else
        edge2vtx_sort = [];
    end
    if ismember('cell2edg', to_check)
        cell2edge_sort = issortedrows(mesh.cell2edg.');
    else
        cell2edge_sort = [];
    end
    sort_all = [vtx_sort, cell2vtx_sort, edge2vtx_sort, cell2edge_sort];
    
    if ~all(sort_all)
        n_fail = length(find(~sort_all));
        check_fail = to_check(~sort_all);
        fprintf(sprintf(['Check failed for:\n', ...
            repmat('- %s\n', 1, n_fail)], check_fail{:}));
    else
        fprintf('All mesh entities are sorted ascending.\n');
    end
end