function [M, TM] = assembleMass(fe, mesh, param, verbosity)
    % Assembles the sparse mass matrix.
    %
    % TODO: implement.
    %
    % a(u,v) = \int_Omega u * \param(x,y) v d(x,y) = ...
    %   \sum_k param_k \abs(\det(B_k)) \sum_l ( w_l ...
    %       \sum_i u_i B_k^(-1) \phi_i({x,y}_l) * ...
    %       \sum_j u_j B_k^(-1) \phi_j({x,y}_l)
    %                              )
    % k   - num simplices
    % l   - num quadrature nodes
    % j,i - num basis functions  
    %
    % Using elemente-wise procedure to set up the global mass matrix.
    %
    % SYNTAX
    %   M = assembleMass(fe, mesh, param[, verbosity])
    %
    % INPUT PARAMETER
    %   fe    ... Struct, including all information to set up Lagrange FE.
    %   mesh  ... Struct, containing mesh information, i.e. coordinates
    %             of vertices and its relation to the triangles and edges.
    %   param      ... Vector of constant cell parameter values.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   M  ... Matrix, representing the mass part of the variational
    %          formulation.
    %   TS ... 3-way-tensor, representing the parameter independent
    %          derivative of M w.r.t. the param vector.
    
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
    [i, j, m, m_TM] = deal(zeros(n_cell * n_entry_loc, 1));
    
    % Set up recurring quantity.
    % Get basis functions for all quadrature nodes referred to
    % the reference simplex.
    quad_eval = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords(:,1), gauss_cords(:,2));      
        
    % Iterate over all simplices.
    for ii = 1:n_cell                   
        % Set up kernel for integral (quadrature summation).
        % By multiplying the vector of basis functions by itselfe
        % using the outer (tensor / dyadic) product the 
        % n_DOF_loc x n_DOF_loc local mass matrix can be obtained in
        % only one step.
        % (note: base_quad_eval is defined as row vector of gradients)
        quad_kern = cellfun(@(x, y) {y * (x.' * x)}, ...
            quad_eval, gauss_weights);
        
        % Apply quadrature summation of the integral over the current
        % simplex.
        % As integral is referred to the reference simplex, the
        % Jacobi-determinat has to be incorporated.
        m_loc = abs(mesh.maps{ii}.detB) * sum(cat(3, quad_kern{:}), 3);
                  
        % Fill up index and value vectors.
        cur_DOF_map = fe.DOF_maps.cell2DOF{ii};
        [i_loc, j_loc] = ndgrid(cur_DOF_map);
        glob_idx_start = ((ii-1) * n_entry_loc) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc - 1;
        i(glob_idx_start:glob_idx_end) = i_loc(:);
        j(glob_idx_start:glob_idx_end) = j_loc(:);
        % Combine constant local cell parameter with the current kernel.
        m(glob_idx_start:glob_idx_end) = param(ii) * m_loc(:);
        % Ignore this parameter value for tensor assembling.
        m_TM(glob_idx_start:glob_idx_end) = m_loc(:);
    end
    
    % Create sparse matrix from index and value vectors.
    % Note, values belonging to the same index pair are automatically
    % summed up by sparse().
    M = sparse(i, j, m, n_DOF_glob, n_DOF_glob);
    
    % Create sparse 3-way-tensor from index and value vectors.
    % Note, that this quantity is independent of the parameter vector. The
    % mass matrix can be obtained by a sparse tensor-times-vector 
    % multiplication:
    %   M = ttv(TM, param, 3);
    size_TM = [n_DOF_glob, n_DOF_glob, fe.sizes.cell];
    index_TM = [i, j, kron((1:fe.sizes.cell).', ones(n_entry_loc, 1))];
    TM = Tensor.Tensor3Coord(size_TM, index_TM, m_TM);
    
    if verbosity
       fprintf('done.\n'); 
    end
end