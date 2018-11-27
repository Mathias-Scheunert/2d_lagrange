function [] = exportObsData2BERT(d_obs, dc_conf, verbosity)
    % Read in DC measurement konfiguration data from BERT formated file.
    
    %% Check input.
    
    assert(ischar(d_obs), ...
        'name - Char, denoting the BERT formatted file, expected.');
    assert(isstruct(dc_conf) &&  isfield(dc_conf, {'map_TX', 'map_RX'}), ...
        'dc_conf - Struct, containing the DC measurement configuration info, expected.'); 
    if nargin < 2
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, dennoting if stats should be print, expected.');
    end

end