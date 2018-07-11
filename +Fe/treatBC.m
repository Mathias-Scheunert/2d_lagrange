function [sol, bnd] = treatBC(fe, mesh, sol, bnd, verbosity)
    % Handles the given BC.
    %
    % SYNTAX
    %   sol = treatBC(fe, mesh, sol, bnd[, verbosity])
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
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct, depending on the BC the system
    %           matrix and or the rhs-vector is appended or reduced,
    %           respectively.
    %   bnd ... Struct, containing the boundary condition information 
    %           appended by bnd DOF information.
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
    
    % Check consistency.
    assert(all(idx_D | idx_N), ...
        'Unknown BC types detected.');
       
    %% Treat Neumann.
    
    % As the entire system might be reduced during the handling of
    % Dirichlet BC, handling of Neumann BC has to processed first.

    if any(idx_N)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_N});
        bnd_N = struct('val', {{bnd.val{idx_N}{~obsolet_BC}}}, ...
                       'DOF', {{bnd.bndDOF.bnd_DOF{~obsolet_BC}}}, ...
                       'param', bnd.param);

        % Process.
        sol = treatNeumann(fe, mesh, sol, bnd_N, verbosity);
    end
    
    %% Treat Dirichlet.

    if any(idx_D)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_D});
        bnd_D = struct('val', vertcat(bnd.val{idx_D}{~obsolet_BC}), ...
                       'DOF', vertcat(bnd.bndDOF.bnd_DOF{~obsolet_BC}));

        % Process.
        sol = treatDirichlet(fe, sol, bnd_D, verbosity);
    end
end

function sol = treatDirichlet(fe, sol, bnd, verbosity)
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
    b_dirichlet = zeros(fe.sizes.DOF ,1);
    if any(bnd.val ~= 0)
        % To couple the inhomogeneouse Dirichlet values with the other DOF
        % the rhs needs to be modyfied.
        % Create a rhs vector which only includes Dirichlet values.
        % Pure edge DOF:
        b_dirichlet(bnd.DOF) = bnd.val;
        
        % Map this vector with the full linear system.
        b_dirichlet_mod = sol.A * b_dirichlet;
        
        % Subtract the result from the original rhs vector
        sol.b = sol.b - b_dirichlet_mod;        
    end
    % Reduce system.
    DOF_req = ~ismember(1:fe.sizes.DOF, bnd.DOF).';
    % RHS vector.
    sol.b = sol.b(DOF_req);
    % System matrix.
    sol.A = sol.A(DOF_req, DOF_req);
    
    %% Append Dirichlet bnd info.
    
    sol.dirichlet.val = b_dirichlet(bnd.DOF);
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
        % Quadrature summation along the edge elements of the boundary.
        % \delta(u) / \delta(n) = g(x, y)
        % \int param(x,y) g(x,y) = \sum_k param_k ...
        %           \sum_j w_j g(x_j,y_j) \sum_i \phi_i(x_j, y_j)
        %
        % k - number of edges (== number of simplices related to those)
        % j - number of quadrature nodes
        % i - number of basis functions

        %% Set up quadrature rule.
        
        % Get quadrature rule for 1D and reshape coordinates such that they
        % can be applied on the basis functions (defined on a 2D reference
        % simplex).
        [gauss_cords, gauss_weights] = Quad.getQuadratureRule(fe.order, 1);
        gauss_cords = [gauss_cords, 0 * gauss_cords].';
        gauss_weights = num2cell(gauss_weights).';
        
        % Set up recurring quantity:
        % Evaluate basis functions for all quadrature nodes referred to
        % the reference simplex at an arbitrary edge (as each edge has got 
        % unit length and all basis functions behave similar w.r.t. the
        % barycentric coordinates it is not required to get the appropriate
        % edge for that).
        basis_eval = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
            gauss_cords(1,:), gauss_cords(2,:)).';
                
        % Iterate over all different Neumann boundaries.
        for kk = 1:length(bnd.val)
        
            %% Set up assembling by obtaining required indices and coords.
            
            % Get only bnd vertices.
            bnd_vtx = bnd.DOF{kk}(bnd.DOF{kk} <= fe.sizes.vtx);
            
            % Get edge index corresponding to given DOFs at Neumann boundary.
            % Note: n_edge + 1 = length(bnd_vtx)
            edge_idx_map = sum(ismember(mesh.edge2vtx, bnd_vtx), 2) > 1;
            edge_idx = find(edge_idx_map);
            n_edge = length(edge_idx);

            % Get indicex of simplex which belong to current edge.
            [~, cell_2_edge_idx_map] = ismember(mesh.cell2edg, edge_idx);
            
            % Get edge coordinates.
            bnd_edge_coo = mesh.edge2cord(edge_idx);

            % Get edge normal vectors.
            % (Reshape output, as only one normal direction can be obtained
            % for bnd edges).
            bnd_edge_n = Mesh.getEdgeNormal(mesh, edge_idx);
            bnd_edge_n = vertcat(bnd_edge_n{:});
            assert(length(bnd_edge_n) == n_edge, ...
                ['Some considered edges are no boundary edge ', ...
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
                fun_hand = RefSol.getConstFunction(bnd.val{kk});
                bnd_val = fun_hand.f;
                gradient = false;
            else
                bnd_val = bnd.val{kk};
                gradient = true;
            end
            
            % Iterate over all edges.
            for ii = 1:n_edge
                                
                % Obtain respective simplex DOFs.
                [cur_cell, ~] = find(cell_2_edge_idx_map == ii);
                cur_DOF_map = fe.DOF_maps.cell2DOF{cur_cell};

                % Get affine maps for current edge (represents a 1D barycentric
                % coordinate system for 2D space).
                Bk = diff(bnd_edge_coo{ii}, 1, 2);
                bk = bnd_edge_coo{ii}(2, :).';
                detBk = sqrt(abs(Bk.' * Bk));

                % Get quadrature nodes position in general form (here gauss 
                % coords in 1D edge related coordinates are required).
                gauss_cords_global = bsxfun(@times, Bk, gauss_cords(1,:)) + bk;               
                
                % Set Neumann BC at the quadrature nodes.
                if gradient
                    % If required, obtain normal derivative.
                    neum_eval = num2cell(arrayfun(@(x, y) ...
                        dot(bnd_val(x, y), bnd_edge_n{ii}), ...
                        gauss_cords_global(1,:), gauss_cords_global(2,:))).';
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