function M = assembleMass(fe, mesh, verbosity)
    % Assembles the sparse mass matrix.
    %
    % continuous:
    % a(u) = \int_Omega |u|^2 d(x,y) = ...
    % Galerkin approx. 
    % (i.e. piece-wise evaluation w.r.t. simplices including the  
    %     numerical quadrature approx for integral evaluation and coord.
    %     shift to reference simplex):
    %   \sum_k
    %       \sum_l ( w_l ...
    %           (1/\abs(\det(B_k) B_k \phi_i({x,y}_l))^T * ...
    %           (1/\abs(\det(B_k) B_k \phi_j({x,y}_l))
    %              )
    %   * \abs(\det(B_k))
    %
    % k   - num simplices
    % l   - num quadrature nodes
    % j,i - num basis functions  
    %
    % Using elemente-wise procedure to set up the global mass matrix.
    %
    % SYNTAX
    %   M = assembleMass(fe, mesh[, verbosity])
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE.
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   M  ... Matrix, representing the mass part of the variational
    %          formulation.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base'})), ...
        'fe - struct, including Lagrange reference element info , expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord', 'maps'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    if nargin < 3
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end

    %% Assemble mass matrix.

    if verbosity
       fprintf('Assemble mass matrix ... '); 
    end
       
    % Get common sizes.
    n_cell = fe.sizes.cell;
    n_DOF_glob = fe.sizes.DOF;
    n_entry_loc = fe.sizes.DOF_loc^2;
    
    % Extract quadrature info from struct.
    gauss_cords = fe.quad.nodes;
    gauss_weights = num2cell(fe.quad.weights);
        
    % Initialize index and value vector for sparse matrix assembling.
    [i, j, m] = deal(zeros(n_cell * n_entry_loc, 1));
        
    % Get basis functions for all quadrature nodes referred to
    % the reference simplex.
    quad_eval = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords(:,1), gauss_cords(:,2));
    
    % Iterate over all simplices.
    for ii = 1:n_cell
        % W.r.t. global coords:
        % To achive normal component of basis functions of adjacent cells
        % to be consistent, check if the related edge normal is parallel 
        % (or antiparallel) oriented with the global normal of the edge
        % (see FeRt.initFiniteElement for its definition).
        %     parallel: Nothing to do. 
        % antiparallel: Switch sign of the basis function.
        % Note: At first, normals are ordered w.r.t. to the mesh.cell2edge 
        % ordering (referred to the relations of cells and edges in global 
        % coords)
        % Get all edges global normal of current cell.
        cur_edge_global_normals = fe.glo_edge_normals(mesh.cell2edg(ii,:), :);
        % Get current edge local normal.
        cur_edge_normals = Mesh.getEdgeNormal(mesh, ...
                               mesh.cell2edg(ii,:), ii);
        cur_edge_normals = cell2mat([cur_edge_normals{:}].');
        % Obtain sign function.
        Phi_sign = dot(cur_edge_global_normals, cur_edge_normals, 2);
        % Change ordering such that sign function can be applied on (local)
        % basis functions for reference simplex.
        Phi_sign = Phi_sign(mesh.loc2glo).';
               
        % Set up fix (related to coordinate transformation) quantity.
        BkTBk = (mesh.maps{ii}.B * 1/abs(mesh.maps{ii}.detB)).' * ...
                (mesh.maps{ii}.B * 1/abs(mesh.maps{ii}.detB));
        
        % Set up kernel for integral (quadrature summation).
        % By multiplying the vector of basis functions by itselfe
        % using the outer (tensor / dyadic) product the 
        % n_DOF_loc x n_DOF_loc local mass matrix can be obtained in
        % only one step.
%         quad_kern = cellfun(@(x, y) {y * (x.' * BkTBk * x)}, ... 
%             quad_eval, gauss_weights);
        quad_kern = cellfun(@(x, y) {...
                        y * ((Phi_sign.*x).' * BkTBk * (Phi_sign.*x))}, ...
                        quad_eval, gauss_weights);
        
        % Apply quadrature summation of the integral over the current
        % simplex.
        % As integral is referred to the reference simplex, the
        % norm of the Jacobi-determinat has to be incorporated.
        m_loc = abs(mesh.maps{ii}.detB) * sum(cat(3, quad_kern{:}), 3);
                  
        % Fill up index and value vectors.
        cur_DOF_map = fe.DOF_maps.cell2DOF{ii};
        [i_loc, j_loc] = ndgrid(cur_DOF_map);
        glob_idx_start = ((ii-1) * n_entry_loc) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc - 1;
        i(glob_idx_start:glob_idx_end) = i_loc(:);
        j(glob_idx_start:glob_idx_end) = j_loc(:);
        m(glob_idx_start:glob_idx_end) = m_loc(:);
    end
    
    % Create sparse matrix from index and value vectors.
    % Note, values belonging to the same index pair are automatically
    % summed up by sparse().
    M = sparse(i, j, m, n_DOF_glob, n_DOF_glob);
    
    if verbosity
       fprintf('done.\n'); 
    end
end