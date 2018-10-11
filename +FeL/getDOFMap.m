function DOF_map = getDOFMap(mesh, fe)
    % Provides mapping of local DOF from reference simplex to a global DOF.
    %
    % SYNTAX
    %   fe = getDOFMap(mesh, fe)
    %
    % INPUT PARAMETER
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %   fe    ... Struct, including affine maps and Lagrange basis
    %             information.
    %
    % OUTPUT PARAMETER
    %   DOF_map ... Struct, containing the DOF map information.
    %
    % REMARKS
    %   The global DOF indices are associated with 
    %   index_of_DOF = [indes_of_verticies, index_of_edges].'
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order'})), ...
        'fe - struct, containing order of Lagrange elements, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, containing cells info, expected.');
    
    %% Define map.
    
    n_vtx = size(mesh.vertices, 1);
    n_cell = size(mesh.cell2vtx, 1);
    switch fe.order
        case 1
            % Set number of DOF.
            n_DOF = n_vtx;
            
            % Collect all DOF for each cell.
            cell2DOF_glo = mesh.cell2vtx;
        case 2
            
            % Set number of DOF.
            n_DOF = n_vtx + size(mesh.edge2vtx, 1);
            
            % Collect all DOF for each cell.
            % [ind_global_vertex-wrt-cell_nodes, ind_global_edge-wrt-cell_edges]
            cell2vert_glo = mesh.cell2vtx;
            cell2edge_glo = mesh.cell2edg; 
            cell2DOF_glo = [cell2vert_glo; n_vtx + cell2edge_glo];
        otherwise
            error('Unsupported order of Lagrangian elements.');
    end
    
    %% Reshape map.
        
    cell2DOF = cell(n_cell, 1);
    switch fe.order
        case 1
            for ii = 1:n_cell
                % Vertex relations.
                glob_vtx = cell2DOF_glo(ii, mesh.loc2glo);
                cell2DOF{ii} = glob_vtx(:);
            end
        case 2
            for ii = 1:n_cell
                % Vertex relations.
                glob_vtx = cell2DOF_glo(ii, mesh.loc2glo);
                % Edge relations.
                %
                % As the edge midpoints are not represented in the original
                % mesh structure, the edge index itself is used for that.
                % As the both, ordering of basis functions on local/reference 
                % simplex and ordering of edges on arbitrary triangle in mesh
                % follows the same principle, both vectors coincide.
                % (See Mesh.appendElementInfo for derivation of the edge 
                % index and definititions in Mesh.getAffineMap.m)
                %                          ind_node, ind_edge
                glob_edg = cell2DOF_glo(n_cell + ii, mesh.loc2glo);
                cell2DOF{ii} = [glob_vtx(:); glob_edg(:)];                
            end
    end

    %% Get DOF coordinates.
    
    switch fe.order
        case 1
            coords = mesh.vertices;
        case 2
            coo_vtx = mesh.vertices;
            coo_edg = cell2mat(cellfun(@(x) {[mean(x(:,1)), mean(x(:,2))]}, ...
                mesh.edge2cord));
            coords = [coo_vtx; coo_edg];            
    end
    
    %% Summarize info.
    
    DOF_map = struct();
    DOF_map.n_DOF = n_DOF;
    DOF_map.cell2DOF = cell2DOF;
    DOF_map.DOF_coo = coords;
end