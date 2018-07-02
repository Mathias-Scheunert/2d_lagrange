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
    % REMARKS
    %   The order of Dirichlet bnd values given in bnd.val (ONLY in case of
    %   length(bnd.val) == 4) is:
    %   bnd.val = [bot; top; left; right]
    %
    % TODO: implement.

    %% Check input
    
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
    if nargin < 5
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    if ~isfield(bnd, 'bndDOF')
        bnd.bndDOF = Fe.getBndDOF(fe, mesh);
    end
    
    %% Treat Dirichlet.
    
    if strcmp(bnd.type, 'dirichlet')
        [sol, bnd] = Fe.treatDirichlet(fe, mesh, sol, bnd, verbosity);
    end
    
    %% Treat Neumann.
    
    if strcmp(bnd.type, 'neumann')
        [sol, ~] = Fe.treatNeumann(fe, mesh, sol, bnd, verbosity);
    end
end