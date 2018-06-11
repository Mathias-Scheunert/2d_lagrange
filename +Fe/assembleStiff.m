function S = assembleStiff(fe, param)
    % Assembles the sparse stiffness matrix.
    %
    % Using elemente-wise procedure to set up the global mass matrix.
    %
    % SYNTAX
    %   S = assembleStiff(fe, param)
    %
    % INPUT PARAMETER
    %   fe    ... Struct, including all information to set up Lagrange FE.
    %   param ... Vector, defining the cell piece-wise constant parameter.
    %
    % OUTPUT PARAMETER
    %   S ... Matrix, representing the stiffness part of the variational
    %         formulation.
    %
    % TODO: implement tensor-based handling (e.g. from toolbox).
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base'})), ...
        'fe - struct, including Lagrange reference element info , expected.');
    assert(length(param) == fe.sizes.cell, ...
        'params - vector with length equal to the number of mesh cells expected.')

    %% Assemble stiffness matrix.
    
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
    [i, j, s] = deal(zeros(n_cell * n_entry_loc, 1));
    
    % Set up recurring quantity.
    % Evaluate gradients of basis functions for all Gauss quadrature nodes 
    % on the reference simplex.
    grad_quad_eval_loc = arrayfun(@(x,y) ...
            {fe.base.grad_Phi(x, y)}, gauss_cords(:,1), gauss_cords(:,2));      
    term2 = cell2mat(grad_quad_eval_loc);
    
    % Iterate over all simplices.
    for ii = 1:n_cell
        % Incorporate inverse of mapping 
        term1 = kron(eye(n_quad_point, n_quad_point), fe.maps{ii}.BinvT);
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
        
        % Combine constant local cell parameter with the integral over the
        % current simplex (quadrature summation).
        % As integral is referred to the reference simplex, the
        % Jacobi-determinat has to be incorporated.
        s_loc = param(ii) * abs(fe.maps{ii}.detB) * ...
            sum(cat(3, quad_kern{:}), 3);
               
        % Fill up index and value vectors.
        cur_DOF_map = fe.DOF_maps.cell2DOF{ii};
        [i_loc, j_loc] = ndgrid(cur_DOF_map);
        glob_idx_start = ((ii-1) * n_entry_loc) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc - 1;
        i(glob_idx_start:glob_idx_end) = i_loc(:);
        j(glob_idx_start:glob_idx_end) = j_loc(:);
        s(glob_idx_start:glob_idx_end) = s_loc(:);
    end
    
    % Create sparse matrix from index and value vectors.
    % Note, values belonging to the same index pair are automatically
    % summed up by sparse().
    S = sparse(i, j, s, n_DOF_glob, n_DOF_glob);
end