function mesh = appendEdgeInfo(mesh)
    % Append element-edge relation.
    % 
    % SYNTAX
    %   mesh = appendEdgeInfo(mesh[, refine])
    %
    % INPUT PARAMETER
    %   mesh   ... Struct, containing mesh information, i.e. coordinates of
    %              vertices and it's relation to the triangles.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, mesh information appended by the edge information.
    %
    % REMARKS
    %  Edge numbering (and orientation) w.r.t. to global vertex numbers:
    %         point2
    %        /      \
    %       /        \
    %      1^         2v
    %     /            \
    %    /_ _ _ 3>_ _ _ \
    % point1          point 3
    %
    % Note that this relation can only be achieved if vertices in cell2vtx
    % are sorted in acending order!

    %% Check input.
    
    assert(isstruct(mesh) && isfield(mesh, 'cell2vtx'), ...
        'mesh - struct, containing vetrex-triangle relation, expected.');

    %% Add edges.
    
    % Make sure that cell2vtx list is orderd ascendingly.
    % (May not be the case if loaded externaly)
    if any(strcmp(mesh.type, {'gmsh_create', 'gmsh_load'}))
        mesh.cell2vtx = sort(mesh.cell2vtx, 2);
    end
    
    % Line number = number of edges.
    % colums 1 and 2 = first and second vertex.
    % (Note: the current order equals the DOF order for second order
    % Lagrange elements by using the Vandermonde matrix for their
    % definition.)
    n_cells = size(mesh.cell2vtx, 1);
    edge_list_full = zeros(3*n_cells, 2);
    for ii = 1:n_cells
        edge_list_full(3*ii-2,:) = [mesh.cell2vtx(ii, 1), mesh.cell2vtx(ii, 2)];
        edge_list_full(3*ii-1,:) = [mesh.cell2vtx(ii, 2), mesh.cell2vtx(ii, 3)];
        edge_list_full(3*ii,:) = [mesh.cell2vtx(ii, 1), mesh.cell2vtx(ii, 3)];
    end
    
    %% Reduce edge list in order to comprise only unique edges.
    
    % Remove multiply occuring edges.
    [edge_list, ~, red2full] = unique(edge_list_full, ...
                                  'rows', 'first');
    
    %% Check consistency with external mesh.
    
    switch mesh.type
        case {'gmsh_create', 'gmsh_load'}
            % Ascendingly sort given bnd edge list.
%             mesh.bnd_edge2vtx = sort(mesh.bnd_edge2vtx, 2);
            sort_bnd_edge2vtx = sort(mesh.bnd_edge2vtx, 2);

            % Find predetermined boundary edges in total edge list;
            [edge2vtx_in_edge_list, idx_edge2vtx_in_edge_list] = ...
                ismember(edge_list, sort_bnd_edge2vtx, 'rows');                
            
            % Check if every boundary edge is included in total edge list.
            assert(length(find(edge2vtx_in_edge_list)) == ...
                size(sort_bnd_edge2vtx, 1), ...
                ['Not all given boundary edges could be found in grid. ', ...
                'Make sure that all (straight) lines are associated ', ...
                'with a physical line within the Gmsh input file.']);
            mesh = rmfield(mesh, 'bnd_edge2vtx');
            
            % Expand mesh struct with mapping from gmsh bnd edge list to
            % total edge list.
            mesh.gmsh_bnd_edge2total_edge = edge2vtx_in_edge_list;
            mesh.gmsh_bnd_edge2total_edge_map = idx_edge2vtx_in_edge_list;
            
        case {'cube', 'rhomb'}
            % Nothing to do.
            
        otherwise
            error('Unknown mesh type.');
    end
    
    %% Summarize information.
    
    mesh.edge2vtx = edge_list;
    mesh.cell2edg = reshape(red2full, 3, n_cells).';
end