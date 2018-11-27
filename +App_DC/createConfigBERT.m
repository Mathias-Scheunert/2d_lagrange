function [dc_conf] = createConfigBERT(ele, type, verbosity)
    % Create DC measurement konfiguration data in BERT format.
    
    %% Check input.
    
    assert(ismatrix(ele) && size(ele, 2) == 2, ...
        'ele - Matrix [n x 2], denoting the electrode positions, expected.');
    assert(ischar(type) && any(strcmp(type, {'', ''})), ...
        'type - Char, denoting the DC measurement configuration, expected.');
    if nargin < 3
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, dennoting if stats should be print, expected.');
    end

    
end


