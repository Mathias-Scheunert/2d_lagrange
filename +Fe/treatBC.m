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
        if verbosity
           fprintf('done.\n'); 
        end
        return;
    end
    
    % Otherwise.
    error('Not supported yet.');
end