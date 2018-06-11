function [sol, bnd] = treatDirichlet(fe, mesh, sol, bnd)
    % Adapts the FE linear system to handle Dirichlet boundary conditions.
    %
    % SYNTAX
    %   fe = treatDirichlet(fe, mesh, sol, bnd)
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
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct (see INPUT PARAMETER).
    %   bnd ... Struct, containing the boundary condition information 
    %           appended by bnd DOF information.
    %
    % REMARKS
    %   The order of Dirichlet bnd values given in bnd.val (ONLY in case of
    %   length(bnd.val) == 4) is:
    %   bnd.val = [bot; top; left; right]
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(sol) && all(isfield(sol, { 'A', 'b'})), ...
        'sol - struct, containing info about linear system to be solved, expected.');
    assert(isstruct(bnd) && all(isfield(bnd, {'type', 'val'})), ...
        'bnd - struct, containing the boundary condition information, expected.');
    assert(strcmp(bnd.type, 'dirichlet'), ...
        'bnd.type - only handling of Dirichlet values are supported here.');
 
    %% Obtain all DOF, belonging to the boundaries.
    
    if isfield(bnd, 'bndFOF')
        bndDOF = bnd.bndDOF;
    else
        bndDOF = Fe.getBndDOF(fe, mesh);
    end
    
    %% Provide consistency for given bnd constraint.
    
    % Check input.
    assert(isvector(bnd.val), ...
        'bnd.val - vector, containing bnd Dirichlet values expected.');
    
    if any(bnd.val ~= 0)
        % Get sizes.
        n_DOF_bot = length(bndDOF.DOF_bot);
        n_DOF_top = length(bndDOF.DOF_top);
        n_DOF_left = length(bndDOF.DOF_left);
        n_DOF_right = length(bndDOF.DOF_right);

        % Handle supported bnd.val formats.
        if length(bnd.val) == 4
            % Uniform Dirichlet constraints per edge.
            % Expected input format:
            %   bnd.val = [bot, top, left, right]
            % Treat DOF of each edge:
            bot = bnd.val(1) + zeros(n_DOF_bot, 1);
            top = bnd.val(2) + zeros(n_DOF_top, 1);
            left = bnd.val(3) + zeros(n_DOF_left, 1);
            right = bnd.val(4) + zeros(n_DOF_right, 1);

            % Treat corner points which will be shared by two domain edges.
            corner_lb = mean(bnd.val([1, 3]));
            corner_rb = mean(bnd.val([1, 4]));
            corner_lt = mean(bnd.val([2, 3]));
            corner_rt = mean(bnd.val([2, 4]));

        elseif length(bnd.val) == length(bndDOF.bnd_DOF)
            % Arbitrary Dirichlet constraints for each DOF per edge.
            bnd_val_split = mat2cell(bnd.val, ...
                [n_DOF_bot, n_DOF_top, n_DOF_left, n_DOF_right], 1);
            bot = bnd_val_split{1}(:);
            top = bnd_val_split{2}(:);
            left = bnd_val_split{3}(:);
            right = bnd_val_split{4}(:);

            % Treat corner points which will be shared by two domain edges.
            warn_corn = false;
            corner_lb = bnd.val(bndDOF.bnd_DOF == bndDOF.vtx_lb);
            if diff(corner_lb) < eps
                corner_lb = corner_lb(1);
            else
                warn_corn = true;
                corner_lb = mean(corner_lb);
            end
            corner_rb = bnd.val(bndDOF.bnd_DOF == bndDOF.vtx_rb);
            if diff(corner_rb) < eps
                corner_rb = corner_rb(1);
            else
                warn_corn = true;
                corner_rb = mean(corner_rb);
            end
            corner_lt = bnd.val(bndDOF.bnd_DOF == bndDOF.vtx_lt);
            if diff(corner_lt) < eps
                corner_lt = corner_lt(1);
            else
                warn_corn = true;
                corner_lt = mean(corner_lt);
            end
            corner_rt = bnd.val(bndDOF.bnd_DOF == bndDOF.vtx_rt);
            if diff(corner_rt) < eps
                corner_rt = corner_rt(1);
            else
                warn_corn = true;
                corner_rt = mean(corner_rt);
            end
            if warn_corn
                waring('Inconsistent Dirichlet values at domain corners are averaged.');
            end
        else
            error(['Unknown Dirichlet value distribution. ', ...
                'Please see Fe.getBndDOF for the boundary DOF composition ', ...
                'and the resulting vector length.']);
        end
    end
    
    %% Reduce the linear system.
       
    % As the already known Dirichlet values and the linear equations for 
    % the respective DOF don't need to be considered in the FE linear 
    % system, they can be excluded.
    b_dirichlet = zeros(fe.sizes.DOF ,1);
    if any(bnd.val ~= 0)
        % To couple the inhomogeneouse Dirichlet values with the other DOF
        % the rhs needs to be modyfied.
        % Create a rhs vector which only includes Dirichlet values.
        % Pure edge DOF:
        b_dirichlet(bndDOF.DOF_bot) = bot;
        b_dirichlet(bndDOF.DOF_top) = top;
        b_dirichlet(bndDOF.DOF_left) = left;
        b_dirichlet(bndDOF.DOF_right) = right;
        % Corner DOF:
        b_dirichlet(bndDOF.vtx_lb) = corner_lb;
        b_dirichlet(bndDOF.vtx_rb) = corner_rb;
        b_dirichlet(bndDOF.vtx_lt) = corner_lt;
        b_dirichlet(bndDOF.vtx_rt) = corner_rt;
        % Map this vector with the full linear system.
        b_dirichlet_mod = sol.A * b_dirichlet;
        % Subtract the result from the original rhs vector
        sol.b = sol.b - b_dirichlet_mod;        
    end
    % Reduce.
    % RHS vector.
    sol.b = sol.b(bndDOF.inner_DOF);
    % System matrix.
    sol.A = sol.A(bndDOF.inner_DOF, bndDOF.inner_DOF);
    
    %% Append Dirichlet bnd info.
    
    bnd.bndDOF = bndDOF;
    sol.dirichlet.val = b_dirichlet(bndDOF.bnd_DOF);
    sol.dirichlet.inner_DOF = bndDOF.inner_DOF;
    sol.dirichlet.bnd_DOF = bndDOF.bnd_DOF;
end