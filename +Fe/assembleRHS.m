function b = assembleRHS(fe, mesh, TX)
    % Assembles the rhs vector for different source types.
    %
    % SYNTAX
    %   b = assembleRHS(fe, mesh, TX)
    %
    % INPUT PARAMETER
    %   fe ... Struct, including all information to set up Lagrange FE.
    %   mesh  ... Struct, containing mesh information, i.e. coordinates
    %             of vertices and its relation to the triangles and edges.
    %   TX ... Struct, containing the source information. 
    %
    % OUTPUT PARAMETER
    %   b ... Vector, representing the given source w.r.t. the DOF.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base', 'maps', 'DOF_maps'})), ...
        'fe - struct, including all information to set up Lagrange FE, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    assert(isstruct(TX) && all(isfield(TX, {'ref_sol'})), ...
        'TX - struct, including source information, expected.');
    
    %% Assemble rhs vector.
    
    % The assembling of the rhs vector follows the same procedure as
    % assembling the mass matrix.
    % I.e. the numerical evaluation of an integral over the single
    % elements/cells is required, whereas the product of the basis
    % functions with a function (describing the source) forms the integral
    % kernel.
    
    % Get common sizes.
    n_cell = fe.sizes.cell;
    n_DOF_glob = fe.sizes.DOF;
    n_entry_loc = fe.sizes.DOF_loc;
    
    % Extract quadrature info from struct.
    gauss_cords = fe.quad.nodes;
    gauss_weights = num2cell(fe.quad.weights.');
        
    % Initialize index and value vector for sparse matrix assembling.
    [i, s] = deal(zeros(n_cell * n_entry_loc, 1));
    
    % Set up recurring quantity.
    % Get basis functions for all quadrature nodes referred to
    % the reference simplex.
    basis_eval = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords(:,1), gauss_cords(:,2)).';
        
    % Iterate over all simplices.
    for ii = 1:n_cell
        
        % Get reference/source function for all quadrature nodes in simplex
        % w.r.t. to the global coordinate system.
        coord_ref = bsxfun(@plus, fe.maps{ii}.B * gauss_cords.', fe.maps{ii}.b);
        fun_eval = arrayfun(@(x, y) {TX.ref_sol.f(x, y)}, ...
                coord_ref(1,:), coord_ref(2,:));
            
        % Set up kernel for integral (quadrature summation).
        % By multiplying the vector of basis functions with the
        % reference/source function.
        quad_kern = cellfun(@(x, y, z) {x * (y * z)}, ...
            basis_eval, fun_eval, gauss_weights);
        
        % Evaluate numerical integration and incorporate Jacobi-determinat 
        % due to mapping back from reference simplex to global coordinates.
        m_loc = abs(fe.maps{ii}.detB) * sum(cat(3, quad_kern{:}), 3);
                  
        % Fill up index and value vectors.
        i_loc = fe.DOF_maps.cell2DOF{ii};
        glob_idx_start = ((ii-1) * n_entry_loc) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc - 1;
        i(glob_idx_start:glob_idx_end) = i_loc(:);
        s(glob_idx_start:glob_idx_end) = m_loc(:);
    end
    
    % Create sparse rhs vector.
    b = sparse(i, 1, s, n_DOF_glob, 1);
end