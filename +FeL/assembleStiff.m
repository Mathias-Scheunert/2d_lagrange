function [S, TS] = assembleStiff(fe, mesh, param, verbosity)
    % Assembles the sparse stiffness matrix.
    %
    % continuous:
    % a(u,v) = \int_Omega \grad(u) * \param(x,y) \grad(v) d(x,y) = ...
    % Galerkin approx. 
    % (i.e. piece-wise evaluation w.r.t. simplices including the  
    %     numerical quadrature approx for integral evaluation and coord.
    %     shift to reference simplex):
    % a(u_i, v_j) = ...
    %   \sum_k param_k ...
    %       \sum_l ( w_l ...
    %           B_k^(-1) \grad(\phi_i({x,y}_l)) * ...
    %           B_k^(-1) \grad(\phi_j({x,y}_l))
    %              )
    %   * \abs(\det(B_k))
    %
    % k   - num simplices
    % l   - num quadrature nodes
    % j,i - num basis functions    
    %
    % Using elemente-wise procedure to set up the global stiffness matrix.
    %
    % SYNTAX
    %   [S, TS] = assembleStiff(fe, mesh, param[, verbosity])
    %
    % INPUT PARAMETER
    %   fe    ... Struct, including all information to set up Lagrange FE.
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %   param ... Vector of constant cell parameter values.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   S  ... Matrix, representing the stiffness part of the
    %          variational formulation.
    %   TS ... 3-way-tensor, representing the parameter independent
    %          derivative of S w.r.t. the param vector.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base'})), ...
        'fe - struct, including Lagrange reference element info , expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord', 'maps'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    assert(isvector(param) && length(param) == size(mesh.cell2vtx, 1), ...
        'param - Vector of constant cell parameter values, expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end

    %% Assemble stiffness matrix.

    if verbosity
       fprintf('Assemble stiffness matrix ... '); 
    end
    
    % Get common sizes.
    n_cell = fe.sizes.cell;
    n_DOF_glob = fe.sizes.DOF;
    n_DOF_loc = fe.sizes.DOF_loc;
    n_quad_point = fe.sizes.quad_point;
    n_entry_loc = n_DOF_loc^2;
    
    % Extract quadrature info from struct.
    gauss_cords = fe.quad.nodes;
    gauss_weights = num2cell(fe.quad.weights);
    
    % Initialize index and value vector for sparse matrix assembling.
    [i, j, s, s_TS] = deal(zeros(n_cell * n_entry_loc, 1));
    
    % Set up recurring quantity.
    % Evaluate gradients of basis functions for all Gauss quadrature nodes 
    % on the reference simplex.
    grad_quad_eval_loc = arrayfun(@(x,y) ...
            {fe.base.grad_Phi(x, y)}, gauss_cords(:,1), gauss_cords(:,2));      
    term2 = cell2mat(grad_quad_eval_loc);
    
    % Iterate over all simplices.
    for ii = 1:n_cell
        % Incorporate inverse of mapping 
        term1 = kron(eye(n_quad_point, n_quad_point), mesh.maps{ii}.BinvT);
        prod = term1 * term2;
        grad_quad_eval = mat2cell(prod, 2 + zeros(n_quad_point, 1), n_DOF_loc);

        % Set up kernel for integral (quadrature summation).
        % By multiplying the vector of basis function gradients by itselfe
        % using the outer (tensor / dyadic) product the 
        % n_DOF_loc x n_DOF_loc local stiffness matrix can be obtained in
        % only one step.
        % (note: base_quad_eval is defined as row vector of gradients)
        quad_kern = cellfun(@(x, y) {y * (x.' * x)}, ...
            grad_quad_eval, gauss_weights);
        
        % Apply quadrature summation of the integral over the current
        % simplex.
        % As integral is referred to the reference simplex, the
        % norm of the Jacobi-determinat has to be incorporated.
        s_loc = abs(mesh.maps{ii}.detB) * sum(cat(3, quad_kern{:}), 3);
               
        % Fill up index and value vectors.
        cur_DOF_map = fe.DOF_maps.cell2DOF{ii};
        [i_loc, j_loc] = ndgrid(cur_DOF_map);
        glob_idx_start = ((ii-1) * n_entry_loc) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc - 1;
        i(glob_idx_start:glob_idx_end) = i_loc(:);
        j(glob_idx_start:glob_idx_end) = j_loc(:);
        % Combine constant local cell parameter with the current kernel.
        s(glob_idx_start:glob_idx_end) = param(ii) * s_loc(:);
        % Ignore this parameter value for tensor assembling.
        s_TS(glob_idx_start:glob_idx_end) = s_loc(:);
    end
    
    % Create sparse matrix from index and value vectors.
    % Note, values belonging to the same index pair are automatically
    % summed up by sparse().
    S = sparse(i, j, s, n_DOF_glob, n_DOF_glob);
    
    % Create sparse 3-way-tensor from index and value vectors.
    % Note, that this quantity is independent of the parameter vector. The
    % stiffness matrix can be obtained by a sparse tensor-times-vector 
    % multiplication:
    %   S = ttv(TS, param, 3);
    size_TS = [n_DOF_glob, n_DOF_glob, fe.sizes.cell];
    index_TS = [i, j, kron((1:fe.sizes.cell).', ones(n_entry_loc, 1))];
    TS = Tensor.Tensor3Coord(size_TS, index_TS, s_TS);
    
    if verbosity
       fprintf('done.\n'); 
    end
end