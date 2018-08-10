function mesh = appendBndInfo(mesh)
    % Append assignment of boundary edges.
    %
    % Do NOT call this function during/after Mesh.refineMeshUniform().
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
    
    % Check if function was called during refinement of gmsh.
    if any(strcmp(mesh.type, {'gmsh_create', 'gmsh_load'}))
        if ~isfield(mesh, 'gmsh_bnd_edge2total_edge') || ...
                isfield(mesh, 'bnd_edge')
            return;
        end
    end
    
    %% Add or handle bnd edge identifier.
    
    switch mesh.type
        case {'cube', 'rhomb'}
            % Get indices from domain boundaries.
            bnd_xmin = cellfun(@(x) all(x(:,1) == mesh.bnd(1)), mesh.edge2cord);
            bnd_xmax = cellfun(@(x) all(x(:,1) == mesh.bnd(2)), mesh.edge2cord);
            bnd_ymin = cellfun(@(x) all(x(:,2) == mesh.bnd(3)), mesh.edge2cord);
            bnd_ymax = cellfun(@(x) all(x(:,2) == mesh.bnd(4)), mesh.edge2cord);
            mesh.bnd_edge = bnd_ymin | bnd_ymax | bnd_xmin | bnd_xmax;
            
        case {'gmsh_create', 'gmsh_load'}
            % Expand gmsh bnd_egde logicals to comprise total edge number.
            mesh.bnd_edge = mesh.gmsh_bnd_edge2total_edge;
            [~, ~, map_order] = find(mesh.gmsh_bnd_edge2total_edge_map);
            
            % Expand domain boundary logicals.
            [bnd_xmin, bnd_xmax, bnd_ymin, bnd_ymax] = ...
                deal(false(size(mesh.edge2vtx, 1), 1));
            bnd_xmin(mesh.gmsh_bnd_edge2total_edge) = mesh.bnd_edge_xmin(map_order);
            bnd_xmax(mesh.gmsh_bnd_edge2total_edge) = mesh.bnd_edge_xmax(map_order);
            bnd_ymin(mesh.gmsh_bnd_edge2total_edge) = mesh.bnd_edge_ymin(map_order);
            bnd_ymax(mesh.gmsh_bnd_edge2total_edge) = mesh.bnd_edge_ymax(map_order);
            
            % Clean up.
            mesh = rmfield(mesh, {'gmsh_bnd_edge2total_edge', ...
                'gmsh_bnd_edge2total_edge_map'});
            
        otherwise
            error('Unkown mesh type');
    end
    
    % Summarize.
    mesh.bnd_edge_xmin = bnd_xmin;
    mesh.bnd_edge_xmax = bnd_xmax;
    mesh.bnd_edge_ymin = bnd_ymin;
    mesh.bnd_edge_ymax = bnd_ymax;
end