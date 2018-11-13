function [fe, sol] = assembleDC2D(mesh, param, fwd_params, verbosity)
    % Set up and assemble the 2D DC linear system and inverse FT info.
    %
    % SYNTAX
    %   [fe, sol] = assembleDC25D(mesh, param, fwd_params, verbosity)
    %
    % INPUT PARAMETER
    %   mesh       ... Struct, containing the mesh information.
    %                  For a detailed description of the content of the
    %                  mesh struct please read header of Mesh.initMesh.
    %   param      ... Vector of constant cell parameter values.
    %   fwd_params ... Struct, containing initial parameters that define
    %                  the forward problem.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   sol     ... Cell [n x 1] of structs for the n wavenumbers.
    %               Each containing the information of the DC problem to be 
    %               solved numerically, i.e. rhs vector, system matrix, 
    %               interpolation operator.
    
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(fwd_params) && all(isfield(fwd_params, ...
        {'TX', 'RX', 'bnd', 'FE_order'})), ...
        ['fwd_params - Struct, containing initial parameters that ', ...
         'define the forward problem., expected.']);
    assert(isvector(param) && length(param) == size(mesh.cell2vtx, 1), ...
        'param - Vector of constant cell parameter values, expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    % Exclude info.
    RX = fwd_params.RX;
    TX = fwd_params.TX;
    bnd = fwd_params.bnd;
    order = fwd_params.FE_order;
    
    %% Set up FE structure.
    
    fe = FeL.initFiniteElement(order, mesh, RX.coo, verbosity);
    bnd = FeL.assignBC(bnd, fe, mesh, param);

    %% Treat 2.5D wavenumber domain handling and set up DC-FE system.

    % Set up invariant rhs vector.
    sol.b = FeL.assembleRHS(fe, mesh, TX, verbosity);

    % Set up invariant system matrix parts.
    A_GradDiv = FeL.assembleStiff(fe, mesh, param, verbosity);
    A_Mass = FeL.assembleMass(fe, mesh, param, verbosity);
    sol.A = A_GradDiv + A_Mass;
    
    % Handle BC..
    sol = FeL.treatBC(fe, mesh, sol, bnd);

    if verbosity
        fprintf('done.\n');
        fprintf('... DC 2D problem initialized.\n \n');
    end
end