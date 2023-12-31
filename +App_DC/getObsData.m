function d_obs = getObsData(u, I, dc_conf)
    % Create observation data from FE solution.
    %
    % SYNTAX
    %   d_obs = getObsData(u, I, dc_conf)
    %
    % INPUT PARAMETER
    %   u       ... Matrix, (block vector) of FE solutions.
    %   I       ... Matrix, interpolation operator which restricts
    %               solutions at DOFs to solution on observation points in
    %               mesh.
    %   dc_conf ... Struct, containing the mapping matrices for superposing
    %               solution for various sources and receiver positions.
    %
    % OUTPUT PARAMETER
    %   d_obs ... Vector, of synthetic data (rhoa at the observations).

    %% Check input.

    assert(ismatrix(u), ...
        'u - Matrix (block solution) from FE forward problem, expected.');
    assert(ismatrix(I), ...
        'I - Matrix of interpolation operator, expected.');
    assert(isstruct(dc_conf) && all(isfield(dc_conf, {'mapTX', 'mapRX'})), ...
        ['dc_conf - Struct, containing the DC measurement ', ...
        'configuration info, expected.']);

    %% Apply mappings.

    % Apply interpolation to map u at DOF to u at all possible RX.
    uRX = I * u;

    % Apply mapping of sources (i.e. superposition of solutions).
    uRX_mapTX = uRX * dc_conf.mapTX;

    % Apply mapping of receivers (i.e. subtracting different solutions).
    u_mapTX_mapRX = diag( dc_conf.mapRX * uRX_mapTX);

    % Transform observed potentials to apparent resistivities by including
    % the configuration factor.
    % Note that solution was obtained for a point source of unit strengh,
    % i.e. normalization w.r.t. current I is not required.
    d_obs = dc_conf.k .* u_mapTX_mapRX;
end

