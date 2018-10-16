function mesh = sortMesh(mesh)
    % Sorts elements of mesh in lexicographic order.
    %
    % The lexicographic / ascending ordering of vertices, cell2vtx, 
    % edge2vtx, cell2edge are required for consistent transition of basis 
    % functions (especially for Raviart-Thomas elements) for two adjacent 
    % cells in mesh.
    %
    % SYNTAX
    %   mesh = sortMesh(mesh)
    %
    % INPUT PARAMETER
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %
    % OUTPUT PARAMETER
    %   mesh ... see INPUT PARAMETER.
    
    %% Check input.
    
    to_sort = {'cell2vtx'};
    opt_to_sort = {'cell2cord'};
    assert(isstruct(mesh) && all(isfield(mesh, to_sort)), ...
        'mesh - Struct, containing mesh info, expected.');
    
    % Add optional field if required.
    to_sort = [to_sort, opt_to_sort(isfield(mesh, opt_to_sort))];
    
    %% Sort elements.
    
    error('Not implemented, yet.');

end