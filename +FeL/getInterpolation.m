function I = getInterpolation(fe, mesh, point)
    % Get map of the DOF to the solution at arbitrary point(s) in the grid.
    %
    % SYNTAX
    %   I = getInterpolation(fe, point)
    %
    % INPUT PARAMETER
    %   fe    ... Struct, including all information to set up Lagrange FE.
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %   point ... Matrix nx2, containing coordinates of observation points.
    %
    % OUTPUT PARAMETER
    %   I ... Matrix, providing the map of the DOF to the solution at point.
    %
    % REMARKS
    %   Usage:
    %       sol = I * u
    %   where
    %       sol - DC potential at point 
    %       u   - solution vector of the FE problem
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base', 'DOF_maps'})), ...
        'fe - struct, including all information to set up Lagrange FE, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord', 'maps'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    assert(isempty(point) || (ismatrix(point) && size(point, 2) == 2), ...
        'point - matrix nx2, containing coordinates, expected.');
    if isempty(point)
        % Leave solution untouched.
        I = eye(fe.sizes.DOF);
        return;
    end
    
    %% Find cell(s) at which the solution is required.
    
    % Get sizes.
    n_point = size(point, 1);    
    
    % Try to identify cells by matlab builtin.
    % Ultra fast variant but may provide NaNs.
    % (May be due to the basic mesh definition or refinemet, inconsistent
    % to delaunay triangulation?)
    cell_idx = tsearchn(mesh.vertices, mesh.cell2vtx, point);
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
        {(Mesh.getAffineMap(x, mesh, point(cell_idx_fail(ii),:)))}, ...
        (1:fe.sizes.cell).');

        % Check if point(s) is/are inside simplex.
        % Note that choosing tol might be quite tricky.
        % Avoid testing points with high number of decimal places.
        tol = pick(3, 0, 1e-5, eps * 1e1);
        cells_fit = cell2mat(cellfun(@(x) {(...
            all(x.xy_ref > -tol, 2) & ...
            all(x.xy_ref <= 1 + tol, 2) & ...
            (sum(x.xy_ref, 2) - 1 < tol)).'}, ...
            maps_fail));
    
        % Obtain cell indices w.r.t to each obervation point.
        % (For multiple hits just take the first cell)
        if isempty(find(cells_fit, 1, 'first'))
            error(['No cell for observation point could be found. ', ...
                'Make sure that observation point lies inside the ', ...
                'modelling domain.']);
        end
        cell_idx(cell_idx_fail(ii)) = find(cells_fit, 1, 'first');
    end
    
    % Get DOF index for respective cells.
    cell_idx = num2cell(cell_idx);
    cells2DOF = cell2mat(cellfun(@(x) {fe.DOF_maps.cell2DOF{x}.'}, cell_idx)).';
    
    % Get affine maps for all found cells w.r.t. the respective points.
    maps = cellfun(@(x, y) {Mesh.getAffineMap(x, mesh, y)}, ...
        cell_idx, mat2cell(point, ones(n_point, 1), 2));
    
    % Get basis function values at observation point(s).
    % (acting as additional scaling in I)
    cur_base = zeros(fe.sizes.DOF_loc, n_point);
    for kk = 1:n_point
        cur_x_ref = maps{kk}.xy_ref(1);
        cur_y_ref = maps{kk}.xy_ref(2);
        cur_base(:, kk) = fe.base.Phi(cur_x_ref, cur_y_ref).';
    end
    
    % Set up interpolation matrix for the linear combination of respective 
    % basis functions.
    n_DOF = fe.sizes.DOF;
    n_DOF_loc = fe.sizes.DOF_loc;
    i = kron((1:n_point).', ones(n_DOF_loc, 1));
    j = cells2DOF(:);
    s = cur_base(:);
    I = sparse(i, j, s, n_point, n_DOF);
end