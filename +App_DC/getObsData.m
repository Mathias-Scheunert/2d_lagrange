function d_obs = getObsData(u, dc_conf)
    % Create observation data from FE solution.

    %% Check input.
    
    assert(ismatrix(u), ...
        'u - Matrix (block solution) from FE forward problem, expected.');
    assert(isstruct(dc_conf) &&  isfield(dc_conf, {'map_TX', 'map_RX'}), ...
        'dc_conf - Struct, containing the DC measurement configuration info, expected.'); 
end

