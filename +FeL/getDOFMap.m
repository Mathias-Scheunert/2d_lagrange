function DOF_map = getDOFMap(mesh, fe)
    % Provides the (re)ordered global DOF list w.r.t. the local2global map.
    %
    % First, all global DOF are listed: DOF = [vertices, n_vtx + edges]
    % Hence, DOF is just a list of all indices:
    %     vertices = 1:n_vtx
    %     edges    = n_vtx + (1:n_edge)
    %
    % Than, for each cell, these indices are shifted in the order that is
    % given by map.loc2glo (see Mesh.getAffineMap).
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

    assert(isstruct(fe) && all(isfield(fe, {'order'})), ...
        'fe - struct, containing order of Lagrange elements, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, containing cells info, expected.');

    %% Collect global DOF from geometry/mesh definition.

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

    %% Reorder global DOF to match loc2glo mapping from affine map definition.

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

    % Note, that the coords have to refer to the initial DOF ordering.
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
