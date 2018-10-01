function [fe, sol, FT_info] = assembleDC25D(mesh, param, fwd_params, verbosity)
    % Set up and assemble the 2.5D DC linear system and inverse FT info.
    %
    % SYNTAX
    %   [fe, sol, FT_info] = assembleDC25D(mesh, param, fwd_params, verbosity)
    %
    % INPUT PARAMETER
    %   mesh       ... Struct, containing mesh information, i.e. 
    %                  coordinatesof vertices and its relation to the
    %                  triangles and edges.
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
    %   FT_info ... Struct, containing the weights and or nodes and the
    %               quadrature type for inverse spatial FT.
    
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(fwd_params) && all(isfield(fwd_params, ...
        {'TX', 'RX', 'bnd', 'FE_order', 'FT_type'})), ...
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
    FT_type = fwd_params.FT_type;
    
    %% Set up FE structure.
    
    fe = FeL.initFiniteElement(order, mesh, RX.coo, verbosity);
    bnd = FeL.assignBC(bnd, fe, mesh, param);

    %% Treat 2.5D wavenumber domain handling and set up DC-FE system.

    % Set up invariant rhs vector.
    rhs = FeL.assembleRHS(fe, mesh, TX, verbosity);

    % Set up invariant system matrix parts.
    A_GradDiv = FeL.assembleStiff(fe, mesh, param, verbosity);
    A_Mass = FeL.assembleMass(fe, mesh, param, verbosity);

    % Get parameter for inverse spatial Fourier transform.
    FT_info = App_DC.getInvFTParam(TX.coo, RX.coo, FT_type);
    
    % Set up FE systems for all wave numbers.
    if verbosity
        fprintf('Assemble linear system and incorporate BC ... '); 
    end
    sol = cell(FT_info.n, 1);
    for ii = 1:(FT_info.n)
        % Set up system matrix
        sol{ii}.A = A_GradDiv + FT_info.k(ii)^2 * A_Mass;

        % Set up rhs.
        sol{ii}.b = (1 / 2) * rhs;

        % Handle boundary conditions.
        sol{ii} = FeL.treatBC(fe, mesh, sol{ii}, bnd);
    end
    if verbosity
        fprintf('done.\n');
        fprintf('... DC 2.5D problem initialized.\n \n');
    end
end