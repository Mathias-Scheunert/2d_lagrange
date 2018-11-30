function dc_conf = getConfigFromBERTData(name, verbosity)
    % Exclude 2D DC configuration info BERT formatted data.
    %
    % The data structure defined in App_DC.readBERTData is assumed without
    % additional checks.
    % Currently only the handling of apparent resistivity (rhoa/Ra), the 
    % resistance (rho/r) or the ratio of voltage (u/U) to current (i/I) is 
    % supported.
    % All data will be stored as apparent resistivity. I.e. rho and U/I
    % will be transformed to that data type. 
    % If rhoa and redundand information (rho, U, I) is given, the latter
    % will be ignored.
    %
    % SYNTAX
    %   dc_conf = getConfigFromBERTData(name, verbosity)
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file name to read from.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   dc_conf ... Struct, containing 
    %               - TX and RX position vectors, 
    %               - mapping matrix for TX (i.e. appropriate summing up 
    %                 the columns in the block rhs for each observation),
    %               - mapping matrix for RX (i.e. appropriate summing up 
    %                 lines in interpolation matrix I for each
    %                 observation),
    %               - configuration factor k,
    %               - observed data vector, defined as apparent resistivity
    %                 rhoa.
    
    %% Check input.
    
    assert(ischar(name), ...
        'name - Char, denoting the BERT formatted file, expected.');
    if nargin < 2
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, dennoting if stats should be print, expected.');
    end

    %% Read out data from BERT formatted file.
    
    all_info = App_DC.readBERTData(name, verbosity);
    
    %% Check data.
    
    % Check if electrode and configuration info is available.
    assert(~isempty(all_info.data_block{1}) && ~isempty(all_info.data_block{2}), ...
        'No electrode of configuration info was found in file.');
    assert(strcmp(all_info.type, '2D'), ...
        'Profiled survey/configuration expected but not found in file.');
    
    %% Prepare routine.
    
    % TODO: may reduce redundancy to App_DC.readBERTData() by removing
    %       alternative tokens there an replace them by standard ones.
    % Predefine standardized header lines.
    head_str = {'Number of electrodes', ...
                'Number of data', ...
                'Number of topo points'};
    n_head = length(head_str);
             
    % Predefine supported standardized data tokens.
    % (These are the subheader column tokens).
    data_token = ...
        {{'x', 'y', 'z'}; ...
         {'a',  'b',  'm',  'n',  'rhoa', 'rho', 'u', 'i'}; ...
         {'x', 'h'} ...
        };
    data_token_alternative = ...
        {{'', '', ''}; ...
         {'c1', 'c2', 'p1', 'p2', 'Ra',   'r',    '',   ''}; ...
         {'', ''} ...
        };
    
    %% Exclude TX/RX info.
    
    % Get/set some parameter.
    n_data = size(all_info.data_block{2}, 1);
    
    % Find columns including index of electrode A and B.
    conf_A = strcmp(data_token{2}{1}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{1}, all_info.line_header{2});
    conf_B = strcmp(data_token{2}{2}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{2}, all_info.line_header{2});
         
    % Find columns including index of electrode A and B.
    conf_M = strcmp(data_token{2}{3}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{3}, all_info.line_header{2});
    conf_N = strcmp(data_token{2}{4}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{4}, all_info.line_header{2});
         
    % Find columns including electrodes x and z coordinates.
    conf_ele_x = strcmp(data_token{1}{1}, all_info.line_header{1}) | ...
                 strcmp(data_token_alternative{1}{1}, all_info.line_header{1});
    conf_ele_z = strcmp(data_token{1}{3}, all_info.line_header{1}) | ...
                 strcmp(data_token_alternative{1}{3}, all_info.line_header{1});
                      
    % Get indices for each TX electrode.
    idx_A = all_info.data_block{2}(:,conf_A);
    idx_B = all_info.data_block{2}(:,conf_B);
    
    % Get indices for each TX electrode.
    idx_M = all_info.data_block{2}(:,conf_M);
    idx_N = all_info.data_block{2}(:,conf_N);
    
    % Determine configuration type.
    % If one of the transmitting electrodes is kept fix, it is assumed that
    % it was used as electrode at infinity, i.e. a pol-dipole measurement
    % is considered (especially if its index == 0).
    % TODO: improve tests. E.g. what if both TX electrodes are kept fix?
    unique_idx_A = unique(idx_A);
    unique_idx_B = unique(idx_B);
    if length(unique_idx_A) == 1 || length(unique_idx_B) == 1
        conf_type = 'pole-dipole';
        if unique_idx_A == 0
            idx_A = [];
        end
        if unique_idx_B == 0
            idx_B = [];
        end
        assert(~isempty([idx_A;  idx_B]), ...
            'All TX electrodes are undefined.');
    else
        assert(length(unique_idx_A) == length(unique_idx_B), ...
            ['Inconsistent electrode indexing obtained. Not sure ', ...
            'which configuration type (pole-dipole or dipole-dipole) ', ...
            'was used.']);
        conf_type = 'dipole-dipole';
    end
        
    % Get positions of all electrodes acting as TX.
    idx_TX = unique([idx_A; idx_B]);
    TX_coo = [all_info.data_block{1}(idx_TX, conf_ele_x), ...
              all_info.data_block{1}(idx_TX, conf_ele_z)];
    
    % Get positions of all electrodes acting as TX.
    idx_RX = unique([idx_M; idx_N]);
    RX_coo = [all_info.data_block{1}(idx_RX, conf_ele_x), ...
              all_info.data_block{1}(idx_RX, conf_ele_z)];
          
    % Set TX strengh (i.e. current I).
    % As only rhoa will be considered as supported data type and U is
    % linear in I, the source strengh can be choosen arbitrarily 
    % (e.g. == 1). 
    TX_val = ones(size(TX_coo, 1), 1);
          
    %% Exclude configuration data.
       
    % Check if supported data tokens can be found.
    data_rhoa = strcmp(data_token{2}{5}, all_info.line_header{2}) | ...
                strcmp(data_token_alternative{2}{5}, all_info.line_header{2});
    data_rho = strcmp(data_token{2}{6}, all_info.line_header{2}) | ...
               strcmp(data_token_alternative{2}{6}, all_info.line_header{2});
    data_U = strcmp(data_token{2}{7}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{7}, all_info.line_header{2});
    data_I = strcmp(data_token{2}{8}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{8}, all_info.line_header{2});
    assert(any([data_rhoa, data_rho, data_U, data_I]), ...
        'No supported data type (rhoa, rho, U, I) was found in file.');

    %% Set up configuration factor k.    
    
    % Obtain RX electrode positions.
    M_coo = [all_info.data_block{1}(idx_M, conf_ele_x), ...
             all_info.data_block{1}(idx_M, conf_ele_z)];
    N_coo = [all_info.data_block{1}(idx_N, conf_ele_x), ...
             all_info.data_block{1}(idx_N, conf_ele_z)];
    getR = @(x) sqrt(x(:,1).^2 + x(:,2).^2);
    switch conf_type
        case 'pole-dipole'
            % Obtain TX electrode positions.
            if isempty(idx_A)
                A_coo = [all_info.data_block{1}(idx_B, conf_ele_x), ...
                         all_info.data_block{1}(idx_B, conf_ele_z)];
            else
                A_coo = [all_info.data_block{1}(idx_A, conf_ele_x), ...
                         all_info.data_block{1}(idx_A, conf_ele_z)];
            end
            
            % Obtain spacings.
            r1 = getR(M_coo - A_coo); % AM
            r2 = getR(N_coo - A_coo); % AN
            
            % Define k.
            k = 1/(2*pi) * 1./(1./r1 - 1./r2);

        case 'dipole-dipole'
            % Obtain TX electrodes positions.
            A_coo = [all_info.data_block{1}(idx_A, conf_ele_x), ...
                     all_info.data_block{1}(idx_A, conf_ele_z)];
            B_coo = [all_info.data_block{1}(idx_B, conf_ele_x), ...
                     all_info.data_block{1}(idx_B, conf_ele_z)];
            
            % Obtain spacings.
            r1 = getR(M_coo - A_coo); % AM
            r2 = getR(N_coo - A_coo); % AN
            r3 = getR(M_coo - B_coo); % BM
            r4 = getR(N_coo - B_coo); % BN
            
            % Define k.
            k = 1/(2*pi) * 1./((1./r1 - 1./r2) - (1./r3 - 1./r4));
    end
    
    %% Exclude vector of observed data (=rhoa per definition).
    
    % Obtain rhoa and check uniqueness.
    if any(data_rhoa)
        rhoa = all_info.data_block{2}(:,data_rhoa);
        
        if any(data_rho)
            % TODO: proof consistency.
            warning('Rhoa and rho provided - ignoring the latter.'); 
        elseif any(data_U) || any(data_I)
            % TODO: proof consistency.
            warning('Rhoa and U or I provided - ignoring the latter.');
        end
    elseif any(data_rho)
        rhoa = all_info.data_block{2}(:,data_rho) .* k;
        
        if any(data_U) || any(data_I)
            % TODO: proof consistency.
            warning('Rho and U or I provided - ignoring the latter.');
        end
    else
        vec_U = all_info.data_block{2}(:,data_U);
        vec_I = all_info.data_block{2}(:,data_I);
        assert(~isempty(length(vec_U)) && (length(vec_U) == length(vec_I)), ...
            'Inconsistency between measured potential and current obtained.');
        
        rhoa = (all_info.data_block{2}(:,data_U) ./ ...
                all_info.data_block{2}(:,data_I)) ...
               .* k;
    end
  
    %% Set up TX mapping.

    %% Set up RX mapping.
    
    %% Summarize.
    
    dc_config = struct();
    dc_config.TX.type = 'point_exact';
    dc_config.TX.coo = TX_coo;
    dc_config.TX.val = TX_val;
    dc_config.RX.coo = RX_coo;
    dc_config.k = [];
    dc_config.rhoa = rhoa;
    dc_config.mapTX = [];
    dc_config.mapRX = [];
end