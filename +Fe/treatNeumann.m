function sol = treatNeumann(fe, mesh, sol, bnd)
    % Adapts the FE linear system to handle Neumann boundary conditions.
    %
    % SYNTAX
    %   fe = treatNeumann(fe, mesh, sol, bnd)
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
    %
    % REMARKS
    %   The order of Neumann bnd values given in bnd.val (ONLY in case of
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
        
    % Check for homogeneous N-BC.
    if all(bnd.val ~= 0)
        % Nothing to do.
        return;
    end
    
    % Otherwise.
    error('Not supported yet.');
    % TODO: implement further.
    
    %% Obtain all DOF, belonging to the boundaries.
    
    if isfield(bnd, 'bndFOF')
        bndDOF = bnd.bndDOF;
    else
        bndDOF = Fe.getBndDOF(fe, mesh);
    end
end