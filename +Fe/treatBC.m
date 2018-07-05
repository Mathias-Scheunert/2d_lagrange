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
            
    %% Treat Dirichlet.

    if any(idx_D)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_D});
        bnd_D = struct('val', vertcat(bnd.val{idx_D}{~obsolet_BC}), ...
                       'DOF', vertcat(bnd.bndDOF.bnd_DOF{~obsolet_BC}));

        % Treat.
        sol = treatDirichlet(fe, sol, bnd_D, verbosity);
    end

    %% Treat Neumann.

    if any(idx_N)
        % Initialize.
        obsolet_BC = cellfun(@isempty, bnd.val{idx_N});
        bnd_N = struct('val', vertcat(bnd.val{idx_N}{~obsolet_BC}), ...
                       'DOF', vertcat(bnd.bndDOF.bnd_DOF{~obsolet_BC}));

        % Treat.
        sol = treatNeumann(fe, mesh, sol, bnd_N, verbosity);
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
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct (see INPUT PARAMETER).
    %
    % TODO: implement inhomogeneouse case.
    
    if verbosity
       fprintf('Incorporate Neumann BC ... '); 
    end
    
    % Check for homogeneous N-BC.
    if all(bnd.val == 0)
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
                
        % Get quadrature rule for 1D and reshape coordinates such that they
        % can be applied on the basis functions (defined on a 2D reference
        % simplex).
        [gauss_cords, gauss_weights] = Quad.getQuadratureRule(fe.order, 1);
        gauss_cords = [gauss_cords, 0 * gauss_cords].';
        gauss_weights = num2cell(gauss_weights).';
        
        % Get only bndDOF indices for the edges where inhomogeneous Neumann
        % BC are set.
        % Note that Fe.assignBC distributes bnd values for every DOF at an
        % bnd (including vertices).
        edge_idx_map = bnd.val & bnd.DOF > fe.sizes.vtx;
        
        % Reduce given bnd quantities and get respective edge index.
        bnd_DOF = bnd.DOF(edge_idx_map);
        edge_idx = bnd_DOF - fe.sizes.vtx;
        n_edge = length(edge_idx);
        
        % Get edge coordinates.
        bnd_edge_coo = mesh.edge2cord(edge_idx);
        
        % Get indices of all cells which belong to those edges.
        [~, cell_2_edge_idx_map] = ismember(mesh.cell2edg, edge_idx);
        cell_2_edge_idx_map = sum(cell_2_edge_idx_map, 2);

        % Initialize index and value vector for sparse vector assembling.
        n_DOF_loc = fe.sizes.DOF_loc;
        [i, s] = deal(zeros(n_edge * n_DOF_loc, 1));

        % Set up recurring quantity:
        % Evaluate basis functions for all quadrature nodes referred to
        % the reference simplex at an arbitrary edge (as each edge has got 
        % unit length and all basis functions behave similar w.r.t. the
        % barycentric coordinates it is not required to get the appropriate
        % edge for that).
        basis_eval = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
            gauss_cords(1,:), gauss_cords(2,:)).';
        
        % Set up bnd val.
        % TODO: fix at Fe.assignBC.m
        bnd_val = arrayfun(@(x) {RefSol.getConstFunction(x)}, ...
            bnd.val(edge_idx_map));
        
        % Iterate over all edges.
        for ii = 1:n_edge

            % Get affine maps for current edge (represents a 1D barycentric
            % coordinate system for 2D space).
            Bk = (bnd_edge_coo{ii}(1, :) - bnd_edge_coo{ii}(2, :)).';
            bk = bnd_edge_coo{ii}(2, :).';
            detBk = sqrt(abs(Bk.' * Bk));
            
            % Get quadrature nodes position in general form.
            gauss_cords_global = (Bk .* gauss_cords) + bk;
            %TODO: fix, here should one obtain 2 different points, I guess
            % mixing up the position of x and y components within the
            % differend variables produces the issue.
            
            % Evaluate Neumann BC at these quadrature nodes.
            neum_eval = num2cell(arrayfun(@(x, y) bnd_val{ii}.f(x, y), ...
                gauss_cords_global(:,1), gauss_cords_global(:,2)));
            
            % Set up integral kernel.
            s = gauss_weights.' * cellfun(@times, neum_eval, basis_eval, ...
                'UniformOutput', false);
            % TODO: fix call: we have 3/6 basis fun_i aand 2 quadnodes_l
            % which have to be summed up here.
            % Store each 3/6 kernels for each edge in s.
            
            % Set up Neumann vector for the linear combination of respective 
            % basis functions.
            i = cells2DOF(:);
            % TODO: fix, get the respective DOF (I guess from bndDOF).
            j = ones(n_point * n_DOF_loc, 1);
        end
        
        % Create sparse Neumann-BC rhs vector.
        n_DOF = fe.sizes.DOF;
        b_N = sparse(i, j, s, n_DOF, 1);
        
        % Summarize original rhs vector and Neumann BC vector.
        sol.b = sol.b + b_N;
    end
    
    if verbosity
       fprintf('done.\n'); 
    end
end