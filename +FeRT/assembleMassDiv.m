function [M, D] = assembleMassDiv(fe, mesh, verbosity)
    % Assembles the sparse mass and divergence matrices.
    %
    % - Mass matrix -
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
    % - Divergence matrix -
    % continuous:
    %  f(v) = \int_Omega \div u d(x,y) = ...
    % Galerkin approx.
    % (i.e. piece-wise evaluation w.r.t. simplices including the
    %     numerical quadrature approx for integral evaluation and coord.
    %     shift to reference simplex):
    %   \sum_k \sum_l ( w_l q_k div(\phi_i({x,y}_l)) 1/\abs(\det(B_k)) )
    %   * \abs(\det(B_k))
    % with
    %   q_k = 1 for all parameters in the triangle k and
    %   q_k = 0 elsewhere
    %
    % k   - num simplices
    % l   - num quadrature nodes
    % j,i - num basis functions
    %
    % Using elemente-wise procedure to set up the global mass matrix.
    %
    % SYNTAX
    %   [M, D] = assembleMassDiv(fe, mesh[, verbosity])
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
    %   D  ... Matrix, representing the divergence part of the
    %          variational formulation.
    %
    % REMARKS
    %   As the basis (q_k) of the parameter space is 1 for all points
    %   within a triangle, i.e. for all quadrature nodes used within
    %   assembling loop, it is not incorporated in the code below.
    % TODO: Implement in order to have a accurate/consistent description.

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

    %% Assemble matrices.

    if verbosity
       fprintf('Assemble mass and divergance matrices ... ');
    end

    % Get common sizes.
    n_cell = fe.sizes.cell;
    n_DOF_glob = fe.sizes.DOF;
    n_entry_loc_M = fe.sizes.DOF_loc^2;
    n_entry_loc_D = fe.sizes.DOF_loc;

    % Extract quadrature info from struct.
    gauss_cords = fe.quad.nodes;
    gauss_weights = num2cell(fe.quad.weights);

    % Initialize index and value vector for sparse matrix assembling.
    [i_M, j_M, m_M] = deal(zeros(n_cell * n_entry_loc_M, 1));
    [j_D, i_D, s_D] = deal(zeros(n_cell * n_entry_loc_D, 1));

    % Set up recurring quantities:
    % Get basis functions for all quadrature nodes referred to
    % the reference simplex and apply gauss weighting.
    quad_eval = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords(:,1), gauss_cords(:,2));
    quad_eval = cellfun(@(x, y) {sqrt(y) * x}, quad_eval, gauss_weights);
    quad_eval = cell2mat(quad_eval);

    % Get divergence of basis functions for all quadrature nodes referred
    % to the reference simplex and apply gauss weighting.
    basis_eval = arrayfun(@(x,y) {fe.base.div_Phi(x, y)}, ...
                    gauss_cords(:,1), gauss_cords(:,2)).';
    basis_eval = cellfun(@(x, y) {y * x}, basis_eval, gauss_weights.');
    basis_eval = cell2mat(basis_eval.');

    % Iterate over all simplices.
    for ii = 1:n_cell
        % W.r.t. global coords:
        % To achive normal component of basis functions of adjacent cells
        % to be consistent, check if the related edge normal is parallel
        % (or antiparallel) oriented with the global normal of the edge
        % (see FeRt.initFiniteElement for its definition).
        %     parallel: Nothing to do.
        % antiparallel: Switch sign of the basis function.
        %
        % As an internal edge is shared by two adjacent cells (each with
        % normals oriented outwards), expoiting the mesh properties:
        % I.e. Consider the cell2edge list and assume the first occuring
        %      cell, sharing an edge to provide the 'global' edge normal.
        %      Check if current cell is the first (= parallel) or second
        %      (=antiparallel) in that list and therefore adjust sign.
        %
        % Get edge indices to search for.
        cur_edge_idx = mesh.cell2edg(ii,:);

        % Search for occurence of current edge(s) in cell2edge list.
        edge_1_occur = mesh.cell2edg == cur_edge_idx(1);
        edge_2_occur = mesh.cell2edg == cur_edge_idx(2);
        edge_3_occur = mesh.cell2edg == cur_edge_idx(3);

        % Obtain indes for first occurence.
        [edge_1_2_cell, ~] = find(edge_1_occur, 1, 'first');
        [edge_2_2_cell, ~] = find(edge_2_occur, 1, 'first');
        [edge_3_2_cell, ~] = find(edge_3_occur, 1, 'first');

        % Derive sign from position w.r.t current cell index.
        Phi_sign = [2 * (edge_1_2_cell == ii) - 1, ...
                    2 * (edge_2_2_cell == ii) - 1, ...
                    2 * (edge_3_2_cell == ii) - 1];

        % As sign function will act on the basis functions in global coords
        % incorporate appropriate DOF mapping.
        Phi_sign = Phi_sign(mesh.loc2glo);

        % Set up fix (related to coordinate transformation) quantity.
        BkTBk = (mesh.maps{ii}.B * 1/abs(mesh.maps{ii}.detB)).' * ...
                (mesh.maps{ii}.B * 1/abs(mesh.maps{ii}.detB));
        BkTBk = kron(diag([1, 1, 1]), BkTBk);
        detBkinv = repmat(1/abs(mesh.maps{ii}.detB), 3, 1).';

        % Set up kernel for integral (quadrature summation) for M.
        % By multiplying the vector of basis functions by itselfe
        % using the outer (tensor / dyadic) product the
        % n_DOF_loc x n_DOF_loc local mass matrix can be obtained in
        % only one step.
        % Evaluate numerical integration and incorporate the norm of the
        % Jacobi-determinat due to mapping back from reference simplex to
        % global coordinates.
        m_loc_M = abs(mesh.maps{ii}.detB) * ...
            ((Phi_sign.*quad_eval).' * BkTBk * (Phi_sign.*quad_eval));

        % Set up kernel for integral (quadrature summation) for D.
        % By multiplying the vector of basis functions with the
        % reference/source function.
        % Evaluate numerical integration and incorporate the norm of the
        % Jacobi-determinat due to mapping back from reference simplex to
        % global coordinates.
        m_loc_D = abs(mesh.maps{ii}.detB) * ...
            sum((Phi_sign.*basis_eval) .* detBkinv);

        % Fill up index and value vectors.
        % w.r.t. M:
        cur_DOF_map = fe.DOF_maps.cell2DOF{ii};
        [i_loc, j_loc] = ndgrid(cur_DOF_map);
        glob_idx_start = ((ii-1) * n_entry_loc_M) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc_M - 1;
        i_M(glob_idx_start:glob_idx_end) = i_loc(:);
        j_M(glob_idx_start:glob_idx_end) = j_loc(:);
        m_M(glob_idx_start:glob_idx_end) = m_loc_M(:);

        % w.r.t. D:
        j_loc = fe.DOF_maps.cell2DOF{ii};
        glob_idx_start = ((ii - 1) * n_entry_loc_D) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc_D - 1;
        j_D(glob_idx_start:glob_idx_end) = j_loc(:);
        i_D(glob_idx_start:glob_idx_end) = ii + (0 * j_loc(:));
        s_D(glob_idx_start:glob_idx_end) = m_loc_D(:);
    end

    % Create sparse matrix from index and value vectors.
    % Note, values belonging to the same index pair are automatically
    % summed up by sparse().
    M = sparse(i_M, j_M, m_M, n_DOF_glob, n_DOF_glob);
    D = sparse(i_D, j_D, s_D, n_cell, n_DOF_glob);

    if verbosity
       fprintf('done.\n');
    end
end
