function [sol, DOF_req] = NtreatBC_dir(fe, mesh, sol, bnd, verbosity)
    % Handles the given BC.
    %
    % SYNTAX
    %   sol = treatBC(fe, mesh, sol, bnd[, verbosity])
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %   sol  ... Struct, containing the information of the current physical
    %            problem to be solved numerically, i.e. rhs vector, system
    %            matrix, interpolation operator.
    %   bnd  ... Struct, containing the boundary condition information.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct, depending on the BC the system
    %           matrix and or the rhs-vector is appended or reduced,
    %           respectively.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.

    %% Check input

    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(sol) && all(isfield(sol, { 'A', 'b'})), ...
        'sol - struct, containing info about linear system to be solved, expected.');
    assert(isstruct(bnd) && all(isfield(bnd, {'type', 'val'})), ...
        'bnd - struct, containing the boundary condition information, expected.');
    assert(iscell(bnd.type), ...
        'bnd.type - cell, containing char(s) of bnd type(s) expected.');
    if nargin < 5
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end

    % Get indices for different BC types.
    idx_D = cellfun(@(x) strcmp(x, 'dirichlet'), bnd.type);
    idx_N = cellfun(@(x) strcmp(x, 'neumann'), bnd.type);
    idx_DtN = cellfun(@(x) strcmp(x, 'dtn'), bnd.type);

    % Check consistency.
    assert(all(idx_D | idx_N | idx_DtN), ...
        'Unknown BC types detected.');

    %% Treat Neumann.

    if any(idx_N)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_N});
        bnd_N = struct('val', {bnd.val{idx_N}(~obsolet_BC).'}, ...
                       'DOF', {bnd.bndDOF.bnd_DOF(~obsolet_BC).'}, ...
                       'param', bnd.param, ...
                       'quad_ord', bnd.quad_ord);

        % Process.
        sol = treatNeumann(fe, mesh, sol, bnd_N, verbosity);
    end

    %% Treat Dirichlet-to-Neumann.

    if any(idx_DtN)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_DtN});
        bnd_DtN = struct('val', {bnd.val{idx_DtN}(~obsolet_BC).'}, ...
                         'DOF', {bnd.bndDOF.bnd_DOF(~obsolet_BC).'}, ...
                         'param', bnd.param, ...
                         'quad_ord', bnd.quad_ord, ...
                         'k', bnd.k);

        % Process.
        sol = treatDtN(fe, mesh, sol, bnd_DtN, verbosity);
    end

    %% Treat Dirichlet.

    % As the entire system might be reduced during the handling of
    % Dirichlet BC, handling has to be processed last.

    if any(idx_D)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_D});
        bnd_D = struct('val', vertcat(bnd.val{idx_D}{~obsolet_BC}), ...
                       'DOF', vertcat(bnd.bndDOF.bnd_DOF{~obsolet_BC}));

        % Process.
        [sol, DOF_req] = treatDirichlet(fe, sol, bnd_D, verbosity);
    end
end

function [sol, DOF_req] = treatDirichlet(fe, sol, bnd, verbosity)
    % Adapts the FE linear system to handle Dirichlet boundary conditions.
    %
    % SYNTAX
    %   fe = treatDirichlet(fe, sol, bnd[, verbosity])
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   sol  ... Struct, containing the information of the current physical
    %            problem to be solved numerically, i.e. rhs vector, system
    %            matrix, interpolation operator.
    %   bnd  ... Struct, containing the boundary condition information.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct (see INPUT PARAMETER).

    %% Reduce the linear system.

    if verbosity
       fprintf('Incorporate Dirichlet BC ... ');
    end

    % As the already known Dirichlet values and the linear equations for
    % the respective DOF don't need to be considered in the FE linear
    % system, they can be excluded.
    b_dirichlet = zeros(size(sol.b));
    if any(bnd.val ~= 0)
        % To couple the inhomogeneouse Dirichlet values with the other DOF
        % the rhs needs to be modyfied.
        % Create rhs vector(s) which only includes Dirichlet values.
        % Pure edge DOF:
        b_dirichlet(bnd.DOF,:) = repmat(bnd.val, 1, size(sol.b,2));

        % Map this vector with the full linear system.
        b_dirichlet_mod = sol.A * b_dirichlet;

        % Subtract the result from the original rhs vector
        sol.b = sol.b - b_dirichlet_mod;
    end
    % Reduce system.
    DOF_req = ~ismember(1:fe.sizes.DOF, bnd.DOF).';
    % RHS vector.
    sol.b = sol.b(DOF_req,:);
    % System matrix.
    sol.A = sol.A(DOF_req, DOF_req);
    % derivative system matrix
    if isfield(sol,'TA')
        %sol.TA_no_bd = sol.TA;
        %sol.TA_index= index_xy(sol.TA, DOF_req, DOF_req);
        sol.TA = Tensor.restrictTensor(sol.TA,DOF_req);
        %sol.TS = Tensor.restrictTensor(sol.TS,DOF_req);
        %sol.TMk = Tensor.restrictTensor(sol.TMk,DOF_req);
    end


    %% Append Dirichlet bnd info.

    sol.dirichlet.val = b_dirichlet(bnd.DOF,1);
    sol.dirichlet.DOF_req = DOF_req;
    sol.dirichlet.bnd_DOF = bnd.DOF;

    if verbosity
       fprintf('done.\n');
    end
end

function sol = treatNeumann(fe, mesh, sol, bnd, verbosity)
    % Adapts the FE linear system to handle Neumann boundary conditions.
    %
    % SYNTAX
    %   fe = treatNeumann(fe, mesh, sol, bnd[, verbosity])
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %   sol  ... Struct, containing the information of the current physical
    %            problem to be solved numerically, i.e. rhs vector, system
    %            matrix, interpolation operator.
    %   bnd  ... Struct, containing the boundary condition information.
    %            In case of a function handle, this handle needs to
    %            evaluate the function gradient!
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct (see INPUT PARAMETER).
    %           Depending of 'type' nothing is changed (homogenous case) or
    %           the rhs vector is expanded (inhomogenoues case).

    if verbosity
       fprintf('Incorporate Neumann BC ... ');
    end

    % Check for homogeneous N-BC.
    if ~any(cellfun(@(x) isa(x, 'function_handle'), bnd.val)) ...
            && all(cellfun(@(x) x == 0, bnd.val))
        % Nothing to do.
    else

        % Assembling is similar to assembling of rhs or interpolation,
        % except the fact that the integrals are only along the edges (1D).
        % See Fe.assembleRHS, Fe.getInterpolation
        %
        % continuous:
        % f(v) = \int_dOmega_N param(x,y) g(x,y) v ds
        % where
        % \delta(u) / \delta(n) = n * g(x,y)
        % Galerkin approx.
        % (i.e. piece-wise evaluation w.r.t. simplices including the
        %     numerical quadrature approx for integral evaluation and coord.
        %     shift to reference simplex):
        % \sum_k param_k ( ...
        %     \sum_j w_j [n_k * g(x_j,y_j)] \phi_i(\hat{x_j, y_j})
        %                ) \det(Bk)_1D
        %
        % where
        %   \det(Bk)_1D = \sqrt(\abs(\det(Bk_1D^T Bk_1D)))
        %
        %   n_k            ... current edge normal vector
        %   g(x_j,y_j)     ... kernel function evaluated at quadrature node
        %   \hat{x_j, y_j} ... quadrature nodes on references coords
        %   {x_j, y_j}     ... quadrature nodes on global coords
        %
        % k - number of edges (== number of simplices related to those)
        % j - number of quadrature nodes
        % i - number of basis functions

        %% Set up quadrature rule.

        % Get quadrature rule for 1D and reshape coordinates such that they
        % can be applied on the basis functions (defined on a 2D reference
        % simplex).
        [gauss_cords, gauss_weights] = Quad.getQuadratureRule(bnd.quad_ord, 1);
        gauss_cords = [gauss_cords, 0 * gauss_cords].';
        gauss_weights = num2cell(gauss_weights).';

        % Set up recurring quantity:
        basis_eval_all = cell(3, 1);
        % Set up gauss nodes for evaluation on reference simplex edges.
        if mesh.loc2glo_orient(1) < 0
            gauss_cords_1 = fliplr([gauss_cords(1,:); gauss_cords(2,:)]);
        else
            gauss_cords_1 = [gauss_cords(1,:); gauss_cords(2,:)];
        end
        if mesh.loc2glo_orient(2) < 0
            gauss_cords_2 = fliplr([1 - gauss_cords(1,:); gauss_cords(1,:)]);
        else
            gauss_cords_2 = [1 - gauss_cords(1,:); gauss_cords(1,:)];
        end
        if mesh.loc2glo_orient(3) < 0
            gauss_cords_3 = rot90(gauss_cords, 2);
        else
            gauss_cords_3 = flipud(gauss_cords);
        end

        % Evaluate all basis functions for all quadrature nodes referred to
        % all edges of the reference simplex (see Mesh.getAffineMap).
        % I.e. basis function evaluation takes place  on the reference
        % coords!
        basis_eval_all{1} = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
            gauss_cords_1(1,:), gauss_cords_1(2,:)).';
        basis_eval_all{2} = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
            gauss_cords_2(1,:), gauss_cords_2(2,:)).';
        basis_eval_all{3} = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
            gauss_cords_3(1,:), gauss_cords_3(2,:)).';

        % Iterate over all different Neumann boundaries.
        for kk = 1:length(bnd.val)

            %% Set up assembling by obtaining required indices and coords.

            % Get only bnd vertices.
            bnd_vtx = bnd.DOF{kk}(bnd.DOF{kk} <= fe.sizes.vtx);

            % Get edge index corresponding to given DOFs at Neumann
            % boundary.
            % Note: n_edge + 1 = length(bnd_vtx)
            edge_idx_map = sum(ismember(mesh.edge2vtx, bnd_vtx), 2) > 1;
            edge_idx = find(edge_idx_map);
            n_edge = length(edge_idx);

            % Get indicex of simplex which belong to current edge.
            [~, cell_2_edge_idx_map] = ismember(mesh.cell2edg, edge_idx);

            % Get edge coordinates.
            bnd_edge_coo = mesh.edge2cord(edge_idx);

            % Get edge normal vectors (global cords).
            bnd_edge_n = Mesh.getEdgeNormal(mesh, edge_idx);
            bnd_edge_n = vertcat(bnd_edge_n{:});
            % Check if more than one normal was obtained
            % (must not be the case for bnd edges).
            assert(length(bnd_edge_n) == n_edge, ...
                ['Some considered edges aren`t boundary edges ', ...
                '(multiple normals obtained).']);

            %% Assemble Neumann rhs vector.

            % Initialize index and value vector for sparse vector assembling.
            n_DOF_loc = fe.sizes.DOF_loc;
            [i, j, s] = deal(zeros(n_edge * n_DOF_loc, 1));
            j = j + 1;  % (as this variable is const. 1 for vector assembing)

            % Set up bnd val:
            if ~isa(bnd.val{kk}, 'function_handle')
                % If not already given as function handle,
                % transform constant values to this shape.
                fun_hand = RefSol.getConst(bnd.val{kk});
                bnd_val = fun_hand.f;
                gradient = false;
            else
                bnd_val = bnd.val{kk};
                gradient = true;
            end

            % Iterate over all edges.
            for ii = 1:n_edge
                % Obtain respective simplex DOFs and current edge index.
                [cur_cell, cur_cell_glob_edge] = ...
                    find(cell_2_edge_idx_map == ii);
                cur_DOF_map = fe.DOF_maps.cell2DOF{cur_cell};

                % Get the local edge index for the current edge.
                cur_cell_loc_edge = mesh.loc2glo == cur_cell_glob_edge;

                % Get the basis function evaluations for the current edge
                % in local coordinates (and local ordering).
                basis_eval = basis_eval_all{cur_cell_loc_edge};

                % Get affine map for current edge (represents a 1D
                % barycentric coordinate system for 2D space).
                % global maps to local
                %  vtx_1    ->   0
                %  vtx_2    ->   1
                Bk = [bnd_edge_coo{ii}(2,1) - bnd_edge_coo{ii}(1,1);
                      bnd_edge_coo{ii}(2,2) - bnd_edge_coo{ii}(1,2)];
                bk = bnd_edge_coo{ii}(1, :).';
                detBk = sqrt(abs(Bk.' * Bk));

                % Get quadrature nodes position in global coords (gauss
                % nodes in 1D edge related coordinates are required).
                gauss_cords_global = bsxfun(@times, Bk, gauss_cords(1,:)) + bk;

                % Set Neumann BC at the quadrature nodes.
                if gradient
                    % Obtain normal derivative of kernel function.
                    % I.e. function evaluation takes place  on the global
                    % coords!
                    % Note:[2x1] or [1x2] vector-shaped output is expected
                    % from bnd_val!
                    try
                    neum_eval = num2cell(arrayfun(@(x, y) ...
                        dot(bnd_val(x, y), bnd_edge_n{ii}), ...
                        gauss_cords_global(1,:), gauss_cords_global(2,:))).';
                    catch ME
                       if strcmp(ME.identifier, 'MATLAB:dot:InputSizeMismatch')
                           msg = ['Check bnd.val set up within the DRIVE '...
                               '- make shure that a function handle '...
                               'for the GRADIENT is provided.'];
                           causeException = MException(...
                               'MATLAB:treatBC:wrongFunctionHandle', msg);
                           ME = addCause(ME, causeException);
                       end
                       rethrow(ME);
                    end
                else
                    neum_eval = num2cell(arrayfun(@(x, y) ...
                        bnd_val(x, y), ...
                        gauss_cords_global(1,:), gauss_cords_global(2,:))).';
                end

                % Set up integral kernel.
                kern = cellfun(@(x, y, z) {(x * y) .* z}, ...
                    gauss_weights, neum_eval, basis_eval);

                % Evaluate quadrature by summing each basis function w.r.t.
                % quadrature nodes and apply Jacobi-determinat for backtrafo
                % of 1D reference coords.
                % Finally incorporate current simplex parameter value.
                kern_eval = bnd.param(cur_cell) * sum(vertcat(kern{:})) .* detBk;

                % Create index vectors.
                glob_idx_start = ((ii-1) * n_DOF_loc) + 1;
                glob_idx_end = glob_idx_start + n_DOF_loc - 1;
                s(glob_idx_start:glob_idx_end) = kern_eval(:);
                i(glob_idx_start:glob_idx_end) = cur_DOF_map(:);
            end

            % Create sparse Neumann-BC rhs vector.
            n_DOF = fe.sizes.DOF;
            b_N = sparse(i, j, s, n_DOF, 1);

            % Summarize original rhs vector and Neumann BC vector.
            sol.b = sol.b + b_N;
        end
    end

    if verbosity
       fprintf('done.\n');
    end
end

function sol = treatDtN(fe, mesh, sol, bnd, verbosity)
    % Adapts the FE linear system to handle D-t-N boundary condition.
    %
    % SYNTAX
    %   fe = treatDtN(fe, mesh, sol, bnd[, verbosity])
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %   sol  ... Struct, containing the information of the current physical
    %            problem to be solved numerically, i.e. rhs vector, system
    %            matrix, interpolation operator.
    %   bnd  ... Struct, containing the boundary condition information.
    %            The contianed handle needs to evaluate the function
    %            gradient!
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct (see INPUT PARAMETER).
    %           Dirichlet-to-Neumann operator adds a contribution to the
    %           system left-hand-side, i.e. the system matrix.

    % Check integration kernel.
    assert(all(cellfun(@(x) isa(x, 'function_handle'), bnd.val)));

    if verbosity
       fprintf('Incorporate Dirichlet-to-Neumann BC ... ');
    end

    % Assembling can be considered as a combination of a Neumann-BC
    % integral including an additional basisfunction, i.e. providing a mass
    % matrix like expression.
    % See FeL.assembleMass.m
    %
    % continuous:
    %   a_R(u,v) = \int_dOmega_N param(x,y) n * g_R(x,y,TX) u v ds
    % with the integral kernel
    %   g_R(x,y,TX) = ((x, y) - TX) (k/r) (K_1(k, r)/K_0(k, r))
    % Galerkin approx.
    % (i.e. piece-wise evaluation w.r.t. simplices including the
    %     numerical quadrature approx for integral evaluation and coord.
    %     shift to reference simplex):
    %   a_R(u_i, v_j) = ...
    %       \sum_k param_k
    %          \sum_l ( w_l [n_k * g_R({x,y}_l, TX, wave_num)]...
    %               \phi_i(\hat{{x,y}}_l) \phi_j(\hat{{x,y}}_l)
    %                  )
    %       \det(Bk)_1D
    %
    % where
    %   \det(Bk)_1D = \sqrt(\abs(\det(Bk_1D^T Bk_1D)))
    %
    %   n_k               ... current edge normal vector
    %   g_R({x,y}_l, ...) ... kernel function evaluated at quadrature node
    %   \hat[{x,y}}_l     ... quadrature nodes on references coords
    %   {x,y}_l           ... quadrature nodes on global coords
    %
    % k   - num simplices
    % l   - num quadrature nodes
    % j,i - num basis functions

    %% Set up quadrature rule.

    % Get quadrature rule for 1D and reshape coordinates such that they
    % can be applied on the basis functions (defined on a 2D reference
    % simplex).
    [gauss_cords, gauss_weights] = Quad.getQuadratureRule(bnd.quad_ord, 1);
    gauss_cords = [gauss_cords, 0 * gauss_cords].';
    gauss_weights = num2cell(gauss_weights).';

    % Set up recurring quantity:
    basis_eval_all = cell(3, 1);
    % Set up gauss nodes for evaluation on reference simplex edges.
    if mesh.loc2glo_orient(1) < 0
        gauss_cords_1 = fliplr([gauss_cords(1,:); gauss_cords(2,:)]);
    else
        gauss_cords_1 = [gauss_cords(1,:); gauss_cords(2,:)];
    end
    if mesh.loc2glo_orient(2) < 0
        gauss_cords_2 = fliplr([1 - gauss_cords(1,:); gauss_cords(1,:)]);
    else
        gauss_cords_2 = [1 - gauss_cords(1,:); gauss_cords(1,:)];
    end
    if mesh.loc2glo_orient(3) < 0
        gauss_cords_3 = rot90(gauss_cords, 2);
    else
        gauss_cords_3 = flipud(gauss_cords);
    end

    % Evaluate all basis functions for all quadrature nodes referred to
    % all edges of the reference simplex (see Mesh.getAffineMap).
    % I.e. basis function evaluation takes place on reference coords!
    basis_eval_all{1} = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords_1(1,:), gauss_cords_1(2,:)).';
    basis_eval_all{2} = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords_2(1,:), gauss_cords_2(2,:)).';
    basis_eval_all{3} = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        gauss_cords_3(1,:), gauss_cords_3(2,:)).';

    % Iterate over all different Neumann boundaries.
    for kk = 1:length(bnd.val)

        %% Set up assembling by obtaining required indices and coords.

        % Get only bnd vertices.
        bnd_vtx = bnd.DOF{kk}(bnd.DOF{kk} <= fe.sizes.vtx);

        % Get edge index corresponding to given DOFs at Neumann
        % boundary.
        % Note: n_edge + 1 = length(bnd_vtx)
        edge_idx_map = sum(ismember(mesh.edge2vtx, bnd_vtx), 2) > 1;
        edge_idx = find(edge_idx_map);
        n_edge = length(edge_idx);

        % Get indicex of simplex which belong to current edge.
        [~, cell_2_edge_idx_map] = ismember(mesh.cell2edg, edge_idx);

        % Get edge coordinates.
        bnd_edge_coo = mesh.edge2cord(edge_idx);

        % Get edge normal vectors (global cords).
        bnd_edge_n = Mesh.getEdgeNormal(mesh, edge_idx);
        bnd_edge_n = vertcat(bnd_edge_n{:});
        % Check if more than one normal was obtained
        % (must not be the case for bnd edges).
        assert(length(bnd_edge_n) == n_edge, ...
            ['Some considered edges aren`t boundary edges ', ...
            '(multiple normals obtained).']);

        %% Assemble Neumann rhs vector.

        % Initialize index and value vector for sparse matrix assembling.
        n_DOF_glob = fe.sizes.DOF;
        n_entry_loc = fe.sizes.DOF_loc^2;
        [i, j, s] = deal(zeros(n_edge * n_entry_loc, 1));

        % Fetch bnd val and set wavenumber.
        bnd_val = bnd.val{kk};
        bnd_val = bnd_val(bnd.k);

        % Iterate over all edges.
        for ii = 1:n_edge
            % Obtain respective simplex DOFs and current edge index.
            [cur_cell, cur_cell_glob_edge] = ...
                find(cell_2_edge_idx_map == ii);
            cur_DOF_map = fe.DOF_maps.cell2DOF{cur_cell};

            % Get the local edge index for the current edge.
            cur_cell_loc_edge = mesh.loc2glo == cur_cell_glob_edge;

            % Get the basis function evaluations for the current edge
            % in local coordinates (and local ordering).
            basis_eval = basis_eval_all{cur_cell_loc_edge};

            % Get affine map for current edge (represents a 1D
            % barycentric coordinate system for 2D space).
            % global maps to local
            %  vtx_1    ->   0
            %  vtx_2    ->   1
            Bk = [bnd_edge_coo{ii}(2,1) - bnd_edge_coo{ii}(1,1);
                  bnd_edge_coo{ii}(2,2) - bnd_edge_coo{ii}(1,2)];
            bk = bnd_edge_coo{ii}(1, :).';
            detBk = sqrt(abs(Bk.' * Bk));

            % Get quadrature nodes position in global coords (gauss
            % nodes in 1D edge related coordinates are required).
            gauss_cords_global = bsxfun(@times, Bk, gauss_cords(1,:)) + bk;

            try
            % Obtain normal derivative of kernel function.
            % I.e. function evaluation takes place  on the global coords!
            % Note:[2x1] or [1x2] vector-shaped output is expected from
            % bnd_val!
            neum_eval = num2cell(arrayfun(@(x, y) ...
                dot(bnd_val(x, y), bnd_edge_n{ii}), ...
                gauss_cords_global(1,:), gauss_cords_global(2,:))).';
            catch ME
               if strcmp(ME.identifier, 'MATLAB:dot:InputSizeMismatch')
                   msg = ['Check bnd.val set up within the DRIVE '...
                       '- make shure that a function handle '...
                       'for the GRADIENT is provided.'];
                   causeException = MException(...
                       'MATLAB:treatBC:wrongFunctionHandle', msg);
                   ME = addCause(ME, causeException);
               end
               rethrow(ME);
            end

            % Set up integral kernel.
            kern = cellfun(@(x, y, z) {(x * y) .* (z.' * z)}, ...
                gauss_weights, neum_eval, basis_eval);

            % Evaluate quadrature by summing each basis function w.r.t.
            % quadrature nodes and apply Jacobi-determinat for backtrafo
            % of 1D reference coords.
            % Finally incorporate current simplex parameter value.
            s_loc = detBk * sum(cat(3, kern{:}), 3);

            % Fill up index and value vectors.
            [i_loc, j_loc] = ndgrid(cur_DOF_map);
            glob_idx_start = ((ii-1) * n_entry_loc) + 1;
            glob_idx_end = glob_idx_start + n_entry_loc - 1;
            i(glob_idx_start:glob_idx_end) = i_loc(:);
            j(glob_idx_start:glob_idx_end) = j_loc(:);
            % Combine constant local cell parameter with the current kernel.
            s(glob_idx_start:glob_idx_end) = bnd.param(cur_cell) * s_loc(:);
        end

        % Create sparse matrix from index and value vectors.
        % Note, values belonging to the same index pair are automatically
        % summed up by sparse().
        R = sparse(i, j, s, n_DOF_glob, n_DOF_glob);

        % Add contribution of current Robin boundary to system matrix.
        sol.A = sol.A + R;
    end

    if verbosity
       fprintf('done.\n');
    end
end
