function mesh = appendBndInfo(mesh)
    % Append assignment of boundary edges.
    % 
    % SYNTAX
    %   mesh = appendBndInfo(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing mesh information, i.e. coordinates of
    %            vertices and its relation to the triangles and edges as 
    %            well as the coordinate-edge relation.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, mesh information appended by the boundary
    %            information.
    %
    % REMARKS
    %   Note, that this procedure will only work for grids which doesn't
    %   contain any cutting areas and which is bordered by a strict
    %   rectangular shape.

    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'edge2cord'})), ...
        'mesh - appended struct expected as input paramerter.');
       
    %% Add bnd edge identifier.
    
    switch mesh.type
        case {'cube', 'rhomb'}
            bnd_bot = cellfun(@(x) all(x(:,2) == mesh.bnd(3)), mesh.edge2cord);
            bnd_top = cellfun(@(x) all(x(:,2) == mesh.bnd(4)), mesh.edge2cord);
            bnd_left = cellfun(@(x) all(x(:,1) == mesh.bnd(1)), mesh.edge2cord);
            bnd_right = cellfun(@(x) all(x(:,1) == mesh.bnd(2)), mesh.edge2cord);
            mesh.bnd_edge = bnd_bot | bnd_top | bnd_left | bnd_right;
            mesh.bnd_edge_bot = bnd_bot;
            mesh.bnd_edge_top = bnd_top;
            mesh.bnd_edge_left = bnd_left;
            mesh.bnd_edge_right = bnd_right;
            
        otherwise
            % Nothing can be done, just check, if required information was
            % provided by the mesh source.
            assert(isfield(mesh, 'bnd'), ...
                ['For external meshes all boundary information need to ', ...
                'be supplied in advance.']');
    end
end
