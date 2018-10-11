function DOF_map = getDOFMap(mesh)
    % Provides mapping of local DOF from reference simplex to a global DOF.
    %
    % SYNTAX
    %   fe = getDOFMap(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OUTPUT PARAMETER
    %   DOF_map ... Struct, containing the DOF map information.
    %
    % REMARKS
    %   The global DOF indices are associated with 
    %   index_of_DOF = [index_of_edges].'
    
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, containing cells info, expected.');
    
    %% Define map.
    
    % Set number of DOF.
    n_DOF = size(mesh.edge2vtx, 1);

    % Collect all DOF for each cell.
    % [ind_global_edge-wrt-cell_edges]
    cell2DOF_glo = mesh.cell2edg;
    
    %% Reshape map.

    n_cell = size(mesh.cell2vtx, 1);
    
    cell2DOF = cell(n_cell, 1);
    for ii = 1:n_cell
        % Edge relations.
        %
        % As the edge midpoints are not represented in the original
        % mesh structure, the edge index itself is used for that.
        % As the both, ordering of basis functions on local/reference 
        % simplex and ordering of edges on arbitrary triangle in mesh
        % follows the same principle, both vectors coincide.
        % (See Mesh.appendElementInfo for derivation of the edge 
        % index and definititions in Mesh.getAffineMap.m)
        %                 ind_node, ind_edge
        glob_edg = cell2DOF_glo(ii, mesh.loc2glo);
        cell2DOF{ii} = glob_edg(:);                
    end

    %% Get DOF coordinates.
    
    coo_edg = cell2mat(cellfun(@(x) {[mean(x(:,1)), mean(x(:,2))]}, ...
        mesh.edge2cord));
    coords = coo_edg;
    
    %% Summarize info.
    
    DOF_map = struct();
    DOF_map.n_DOF = n_DOF;
    DOF_map.cell2DOF = cell2DOF;
    DOF_map.DOF_coo = coords;
end