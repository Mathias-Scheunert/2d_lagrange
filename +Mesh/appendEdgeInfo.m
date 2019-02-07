function mesh = appendEdgeInfo(mesh)
    % Append element-edge relation.
    % 
    % SYNTAX
    %   mesh = appendEdgeInfo(mesh[, refine])
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            -> no edge information contained.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
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
    %  The same relation is used on local coords (see Mesh.getAffineMap)
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
    n_cells = size(mesh.cell2vtx, 1);
    edge_list_full = zeros(3*n_cells, 2);
    for ii = 1:n_cells
        edge_list_full(3*ii-2,:) = [mesh.cell2vtx(ii, 1), mesh.cell2vtx(ii, 2)];
        edge_list_full(3*ii-1,:) = [mesh.cell2vtx(ii, 2), mesh.cell2vtx(ii, 3)];
        edge_list_full(3*ii,:)   = [mesh.cell2vtx(ii, 1), mesh.cell2vtx(ii, 3)];
    end
    
    %% Reduce edge list in order to comprise only unique edges.
    
    % Remove multiply occuring edges.
    [edge_list, ~, red2full] = unique(edge_list_full, ...
                                  'rows', 'first');
    
    %% Handle external mesh.
    
    switch mesh.type
        case {'gmsh_create', 'gmsh_load'}
            % As Gmsh only provides those edges, which are part of
            % predefined boundary segments, they need to be related to the
            % total list of edges in mesh.
            
            % Ascendingly sort given bnd edge list.
            sort_bnd_edge2vtx = sort(mesh.bnd_edge2vtx, 2);
            check_gmsh_edges = sort_bnd_edge2vtx(:,1) == sort_bnd_edge2vtx(:,2);
            if ~isempty(find(check_gmsh_edges, 1))
                corrupt_edges = find(check_gmsh_edges);
                warning(sprintf(['Corrupt edge detected: edge ', ...
                    repmat('%d ', 1, length(corrupt_edges)), ...
                    'has zero length. Ignoring that edge.'], ...
                    corrupt_edges));
                mesh.bnd_edge2vtx(corrupt_edges, :) = [];
                sort_bnd_edge2vtx(corrupt_edges, :) = [];
            end

            % Find predetermined boundary edges in total edge list;
            [bnd_edge2vtx_in_edge_list, idx_bnd_edge2vtx_in_edge_list] = ...
                ismember(edge_list, sort_bnd_edge2vtx, 'rows');                
            
            % Check if every boundary edge is included in total edge list.
            % TODO: can we just remove multiples here?
            sort_bnd_edge2vtx = unique(sort_bnd_edge2vtx, 'rows');
            assert(length(find(bnd_edge2vtx_in_edge_list)) == ...
                size(sort_bnd_edge2vtx, 1), ...
                ['Not all given boundary edges could be found in grid. ', ...
                'Make sure that all (straight) lines are associated ', ...
                'with a physical line within the Gmsh input file.']);
            mesh = rmfield(mesh, 'bnd_edge2vtx');
            
            % Expand gmsh bnd_edge logicals to comprise total edge number.
            mesh.bnd_edge = bnd_edge2vtx_in_edge_list;
            [~, ~, map_order] = find(idx_bnd_edge2vtx_in_edge_list);
                        
            % Expand domain boundary ids.
            bnd_edge_part = zeros(length(mesh.bnd_edge), 1);
            bnd_edge_part(bnd_edge2vtx_in_edge_list) = ...
                mesh.bnd_edge_part(map_order);
            mesh.bnd_edge_part = bnd_edge_part;
            
        case {'cube'}
            % Nothing to do.
            
        otherwise
            error('Unknown mesh type.');
    end
    
    %% Summarize information.
    
    mesh.edge2vtx = edge_list;
    mesh.cell2edg = reshape(red2full, 3, n_cells).';
end