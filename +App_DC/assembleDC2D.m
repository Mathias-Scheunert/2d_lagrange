function [fe, sol] = assembleDC2D(mesh, param, fwd_params, verbosity)
    % Set up and assemble the 2D DC linear system and inverse FT info.
    %
    % Problem in 2D:
    %      x = [x, y]
    %   \phi = \phi(x)
    %           
    %   -\div(\sigma\grad(\phi)) = I \dirac(x_0) in Omega
    %                       \phi = phi_1         at d_Omega_1 (left)
    %                       \phi = phi_2         at d_Omega_2 (right)
    %               d_\phi / d_n = 0             at d_Omega_3 (top, bottom)
    %
    % 2D Variational problem:
    %   a(u,v) = \int_Omega \grad(\phi') * \sigma \grad(v) + ...
    %                \int_Omega \phi' * \sigma v
    %   f(v)   = I \int_{Omega} \dirac(x_0) v
    %
    % SYNTAX
    %   [fe, sol] = assembleDC2D(mesh, param, fwd_params, verbosity)
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
    %   fe    ... Struct, including all information to set up Lagrange FE.
    %   sol   ... Cell [n x 1] of structs for the n wavenumbers.
    %             Each containing the information of the DC problem to be 
    %             solved numerically, i.e. rhs vector, system matrix, 
    %             interpolation operator.
    
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

    %% Set up 2D DC-FE system.

    % Set up invariant rhs vector.
    sol.b = FeL.assembleRHS(fe, mesh, TX, verbosity);

    % Set up invariant system matrix parts.
    sol.A = FeL.assembleStiff(fe, mesh, param, verbosity);
    
    % Handle BC.
    sol = FeL.treatBC(fe, mesh, sol, bnd);

    if verbosity
        fprintf('done.\n');
        fprintf('... DC 2D problem initialized.\n \n');
    end
end