function u = solveFwd(sol, fe, verbosity)
    % Provides the solution of the FE linear system.
    %
    % SYNTAX
    %   u = solveFwd(fe, verbosity)  
    %
    % INPUT PARAMETER
    %   sol  ... Struct, containing the information of the current physical
    %            problem to be solved numerically, i.e. rhs vector, system 
    %            matrix, interpolation operator.
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %
    % OUTPUT PARAMETER
    %   u         ... Vector, denoting the solution of the FE problem, i.e.
    %                 the coefficients for the linear combination of the FE
    %                 basis functions.
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    
    %% Check input.
    
    assert(isstruct(sol) && all(isfield(sol, { 'A', 'b'})), ...
        'sol - struct, containing info about linear system to be solved, expected.');
    assert(isstruct(fe) && all(isfield(fe, {'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    if nargin < 3
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    %% Obtain solution of fwd problem.
    
    if isfield(sol, 'dirichlet')
        % Obtain solution for reduced system.
        if verbosity
           fprintf('Solve via ''\\'' ... ');
        end
        u_red = sol.A \ sol.b;
        if verbosity
           fprintf('done.\n'); 
        end
        
        % Construct full solution vector (containing the Dirichlet values).
        if verbosity
           fprintf('Treat Dirichlet DOFs ... '); 
        end
        u = zeros(fe.sizes.DOF, size(sol.b, 2));
        u(sol.dirichlet.DOF_req,:) = u_red;
        u(sol.dirichlet.bnd_DOF,:) = repmat(sol.dirichlet.val, 1, size(sol.b, 2));
        if verbosity
           fprintf('done.\n');
        end
        
    else
        % Obtain solution for the full system.
        if verbosity
           fprintf('Solve via ''\\'' ... ');
        end
        u = full(sol.A \ sol.b);
        if verbosity
           fprintf('done.\n'); 
        end
    end

    if verbosity
        fprintf('... FWP solved.\n \n');
    end
end