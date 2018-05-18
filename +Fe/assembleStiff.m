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
    
    % Extract quadrature info from struct.
    gauss_cords = fe.quad.nodes;
    gauss_weights = fe.quad.weights;
    
    % Initialize global stiffness matrix.
    S = zeros(n_DOF_glob, n_DOF_glob);
    
    % Set up recurring quantity.
    grad_quad_eval_loc = arrayfun(@(x) ...
            {fe.base.grad_Phi(gauss_cords(x,1), gauss_cords(x,2))}, ...
            (1:n_quad_point).');        
    term2 = cell2mat(grad_quad_eval_loc);
    
    % Iterate over all simplices.
    for ii = 1:n_cell
        % Get basis function gradients for all quadrature nodes referred to
        % the reference simplex by incorporating the inverse mapping.
%         % Slow variant.
%         grad_quad_eval = cellfun(@(x) {fe.maps{ii}.BinvT * x},  grad_quad_eval_loc);
        % Fast variant.
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
            grad_quad_eval, num2cell(gauss_weights));
        
        % Combine constant local cell parameter with the integral over the
        % current simplex (quadrature summation).
        % As integral is referred to the reference simplex, the
        % Jacobi-determinat has to be incorporated.
        A_loc = param(ii) * abs(fe.maps{ii}.detB) * ...
            sum(cat(3, quad_kern{:}), 3);
        
        % Map entries of local to entries of the global stiffness matrix.
        % Therefore, local entries referring to the same global entry are
        % summed up.
        cur_DOF_map = fe.DOF_maps.cell2DOF{ii};
        S(cur_DOF_map, cur_DOF_map) = S(cur_DOF_map, cur_DOF_map) + A_loc;
    end
    
    % Convert to sparse matrix.
    % TODO: implement tensor-based handling (e.g. from toolbox).
    S = sparse(S);
end