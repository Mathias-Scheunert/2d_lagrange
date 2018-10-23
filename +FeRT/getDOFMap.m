function DOF_map = getDOFMap(mesh)
    % Provides the (re)ordered global DOF list w.r.t. the local2global map.
    %
    % First, all global DOF are listed: DOF = [edges]
    % Hence, DOF is just a list of all indices:
    %     edges = (1:n_edge)
    %
    % Than, for each cell, these indices are shifted in the order that is 
    % given by map.loc2glo (see Mesh.getAffineMap).
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
    %   The mapping relation will be used within the assembling routines:
    %
    %       local DOF, i.e. the basis functions in local coords
    %       related to local positions, i.e. to vertices or edge midpoints
    %       in the reference simplex
    %                                |
    %                                v
    %                     vertices        edges
    %                    [1, 2, 3,       1, 2, 3]
    %
    %                             map to
    %                                |
    %                                v
    %       global DOF, i.e. the basis function in global coords
    %       related to global positions, i.e. to vertices or edge midpoints
    %       in the mesh
    %                     vertices        edges
    %                 [[mesh.loc2glo], [mesh.loc2glo]]
    
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, containing cells info, expected.');
    
    %% Collect global DOF from geometry/mesh definition.
    
    % Set number of DOF.
    n_DOF = size(mesh.edge2vtx, 1);

    % Collect all DOF for each cell.
    % [ind_global_edge-wrt-cell_edges]
    cell2DOF_glo = mesh.cell2edg;
    
    %% Reorder global DOF to match loc2glo mapping from affine map definition.

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
    
    % Note, that the coords have to refer to the initial DOF ordering.
    coo_edg = cell2mat(cellfun(@(x) {[mean(x(:,1)), mean(x(:,2))]}, ...
        mesh.edge2cord));
    coords = coo_edg;
    
    %% Summarize info.
    
    DOF_map = struct();
    DOF_map.n_DOF = n_DOF;
    DOF_map.cell2DOF = cell2DOF;
    DOF_map.DOF_coo = coords;
end