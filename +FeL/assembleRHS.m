function b = assembleRHS(fe, mesh, TX, verbosity)
    % Assembles the rhs vector for different source types.
    %
    % continuous:
    % f(v) = \int_Omega f v d(x,y) = ...
    % Galerkin approx. 
    % (i.e. piece-wise evaluation w.r.t. simplices including the  
    %     numerical quadrature approx for integral evaluation and coord.
    %     shift to reference simplex):
    % f(v_i) = ...
    %   \sum_k \sum_l ( w_l ( ...
    %              \phi_i({x,y}_l) f({x,y}_l))
    %                 )
    %   * \abs(\det(B_k))
    %
    % k - num simplices
    % l - num quadrature nodes
    % i - num basis functions   
    %
    % SYNTAX
    %   b = assembleRHS(fe, mesh, TX[, verbosity])
    %
    % INPUT PARAMETER
    %   fe    ... Struct, including all information to set up Lagrange FE.
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %   TX    ... Struct, containing the source information.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   b ... Vector, representing the given source w.r.t. the DOF.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base', 'DOF_maps'})), ...
        'fe - struct, including all information to set up Lagrange FE, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord', 'maps'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    assert(isstruct(TX), ...
        'TX - struct, including source information, expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    %% Assemble rhs vector.

    if verbosity
       fprintf('Assemble rhs ... '); 
    end
    
    switch TX.type
        case {'reference', 'point_approx'}
            assert(all(isfield(TX, {'ref_sol'})), ...
                'TX.ref_sol - function handle to ref. solution, expected.');
            getRHS = @getFunctionRHS;
            
        case 'point_exact'
            getRHS = @getDistributionRHS;
            
        otherwise
            error('Unknown soure type.');
    end
    
    b = getRHS(fe, mesh, TX);
    
    if verbosity
       fprintf('done.\n'); 
    end
end

function b = getDistributionRHS(fe, mesh, TX)
    % The assembling of the rhs vector for a distribution rhs (Dirac) 
    % follows the same procedure as assembling the interpolation operator.
    % I.e. the evaluation of basisfunctions at the source location whereas
    % the source strengh acts as scaling.

    % Get common sizes.
    n_point = size(TX.coo, 1);
    
    % Try to identify cells belonging to TX.pos by matlab builtin.
    cell_idx = tsearchn(mesh.vertices, mesh.cell2vtx, TX.coo);
    cell_idx_fail = isnan(cell_idx);
    cell_idx_fail = find(cell_idx_fail);
    n_cell_idx_fail = length(cell_idx_fail);
    if n_cell_idx_fail ~= 0
        warning('tsearchn failed to identify cells for some RX points.');
    end
   
    for ii = 1:n_cell_idx_fail
        % Use 'own' functionality to obtain cell index for corrupt points.
        % Get point coordinates w.r.t. the reference simplex.
        maps_fail = arrayfun(@(x) ...
        {(Mesh.getAffineMap(x, mesh, TX(cell_idx_fail(ii),:)))}, ...
        (1:fe.sizes.cell).');

        % Check if point(s) is/are inside simplex.
        tol = pick(2, 0, eps * 1e1);
        cells_fit = cell2mat(cellfun(@(x) {(...
            all(x.xy_ref > -tol, 2) & ...
            all(x.xy_ref <= 1 + tol, 2) & ...
            (sum(x.xy_ref, 2) - 1 < tol)).'}, ...
            maps_fail));
    
        % Obtain cell indices w.r.t to each TX point.
        % (For multiple hits just take the first cell)
        cell_idx(cell_idx_fail(ii)) = find(cells_fit, 1, 'first');
        if isempty(cell_idx(cell_idx_fail(ii)))
            error('No cell for observation point could be found.');
        end
    end
    
    % Get DOF index for respective cells.
    cell_idx = num2cell(cell_idx);
    cells2DOF = cell2mat(cellfun(@(x) {fe.DOF_maps.cell2DOF{x}.'}, cell_idx)).';
    
    % Get affine maps for all found cells w.r.t. the respective points.
    maps = cellfun(@(x, y) {Mesh.getAffineMap(x, mesh, y)}, ...
        cell_idx, mat2cell(TX.coo, ones(n_point, 1), 2));
    
    % Evaluate basis functions at TX point(s) and scale it.
    cur_base = zeros(fe.sizes.DOF_loc, n_point);
    for kk = 1:n_point
        cur_x_ref = maps{kk}.xy_ref(1);
        cur_y_ref = maps{kk}.xy_ref(2);
        cur_base(:, kk) = TX.val(kk) * fe.base.Phi(cur_x_ref, cur_y_ref).';
    end
    
    % Set up rhs vector for the linear combination of respective 
    % basis functions.
    n_DOF = fe.sizes.DOF;
    n_DOF_loc = fe.sizes.DOF_loc;
    i = cells2DOF(:);
    j = ones(n_point * n_DOF_loc, 1);
    s = cur_base(:);
    b = sparse(i, j, s, n_DOF, 1);
end

function b = getFunctionRHS(fe, mesh, TX)
    % The assembling of the rhs vector for a function-like rhs follows the 
    % same procedure as assembling the mass matrix.
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
        coord_ref = bsxfun(@plus, mesh.maps{ii}.B * gauss_cords.', mesh.maps{ii}.b);
        fun_eval = arrayfun(@(x, y) {TX.ref_sol.f(x, y)}, ...
                coord_ref(1,:), coord_ref(2,:));
            
        % Set up kernel for integral (quadrature summation).
        % By multiplying the vector of basis functions with the
        % reference/source function.
        quad_kern = cellfun(@(x, y, z) {x * (y * z)}, ...
            basis_eval, fun_eval, gauss_weights);
        
        % Evaluate numerical integration and incorporate the norm of the 
        % Jacobi-determinat due to mapping back from reference simplex to 
        % global coordinates.
        m_loc = abs(mesh.maps{ii}.detB) * sum(cat(3, quad_kern{:}), 3);
                  
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