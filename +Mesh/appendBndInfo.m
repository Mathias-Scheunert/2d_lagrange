function mesh = appendBndInfo(mesh)
    % Append assignment of boundary edges.
    %
    % Do NOT call this function during/after Mesh.refineMeshUniform().
    % 
    % SYNTAX
    %   mesh = appendBndInfo(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            -> no edge2boundary information contained.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, mesh information appended by the boundary
    %            information.

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
            
            % Summarize.
            mesh.bnd_edge_part_name = {'xmin', 'xmax', 'ymin', 'ymax'};
            mesh.bnd_edge_part = [bnd_xmin, bnd_xmax, bnd_ymin, bnd_ymax] * ...
                                   (1:length(mesh.bnd_edge_part_name)).';
            
        case {'gmsh_create', 'gmsh_load'}
            % Check input.
            assert(isfield(mesh, 'gmsh_bnd_edge2total_edge'), ...
                ['Mapping bewteen the boundary edges given by Gmsh ', ...
                'and the entire edge list is missing. Make sure to ', ...
                'apply Mesh.appendEdgeInfo() after loading a Gmsh mesh.']);
            
            % Expand gmsh bnd_edge logicals to comprise total edge number.
            mesh.bnd_edge = mesh.gmsh_bnd_edge2total_edge;
            [~, ~, map_order] = find(mesh.gmsh_bnd_edge2total_edge_map);
                        
            % Expand domain boundary ids.
            bnd_edge_part = zeros(length(mesh.bnd_edge), 1);
            bnd_edge_part(mesh.gmsh_bnd_edge2total_edge) = ...
                mesh.bnd_edge_part(map_order);
            mesh.bnd_edge_part = bnd_edge_part;

            % Clean up.
            mesh = rmfield(mesh, {'gmsh_bnd_edge2total_edge', ...
                                  'gmsh_bnd_edge2total_edge_map'});
            
        otherwise
            error('Unkown mesh type');
    end
end