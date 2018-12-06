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
                 
    % Predefine supported standardized data tokens.
    % (These are the subheader column tokens).
    data_token = ...
        {{'x', 'y', 'z'}; ...
         {'a',  'b',  'm',  'n',  'rhoa', 'rho', 'u', 'i', 'k'}; ...
         {'x', 'h'} ...
        };
    data_token_alternative = ...
        {{'', '', ''}; ...
         {'c1', 'c2', 'p1', 'p2', 'Ra',   'r',   '',  '' , ''}; ...
         {'', ''} ...
        };
    
    %% Exclude TX/RX info.
    
    if verbosity
       fprintf('Generate measurement configuration info from BERT file ... '); 
    end
    % Get/set some parameter.
    n_data = size(all_info.data_block{2}, 1);
    
    % Find columns including index of electrode A and B.
    conf_A = strcmpi(data_token{2}{1}, all_info.line_header{2}) | ...
             strcmpi(data_token_alternative{2}{1}, all_info.line_header{2});
    conf_B = strcmpi(data_token{2}{2}, all_info.line_header{2}) | ...
             strcmpi(data_token_alternative{2}{2}, all_info.line_header{2});
         
    % Find columns including index of electrode A and B.
    conf_M = strcmpi(data_token{2}{3}, all_info.line_header{2}) | ...
             strcmpi(data_token_alternative{2}{3}, all_info.line_header{2});
    conf_N = strcmpi(data_token{2}{4}, all_info.line_header{2}) | ...
             strcmpi(data_token_alternative{2}{4}, all_info.line_header{2});
         
    % Find columns including electrodes x and z coordinates.
    conf_ele_x = strcmpi(data_token{1}{1}, all_info.line_header{1}) | ...
                 strcmpi(data_token_alternative{1}{1}, all_info.line_header{1});
    conf_ele_y = strcmpi(data_token{1}{2}, all_info.line_header{1}) | ...
                 strcmpi(data_token_alternative{1}{2}, all_info.line_header{1});
    conf_ele_z = strcmpi(data_token{1}{3}, all_info.line_header{1}) | ...
                 strcmpi(data_token_alternative{1}{3}, all_info.line_header{1});
    if any(conf_ele_y) && ~all(all_info.data_block{1}(:,conf_ele_y) == 0)
        warning(['Data set comprises data from 3D coordinate system. ', ...
            'Only 2D problem supported yet. Ingnoring y-coordinate.']);
    end
                      
    % Get indices for each TX electrode.
    idx_A = all_info.data_block{2}(:,conf_A);
    idx_B = all_info.data_block{2}(:,conf_B);
    
    % Get indices for each TX electrode.
    idx_M = all_info.data_block{2}(:,conf_M);
    idx_N = all_info.data_block{2}(:,conf_N);
    
    % Determine configuration type.
    % If one of the transmitting electrodes have value 0, it is expected
    % to be at infinity, i.e. a pol-dipole measurement is considered.
    idx_A_at_inf = idx_A == 0;
    idx_B_at_inf = idx_B == 0;
        
    % Get positions of all electrodes acting as TX.
    % (Exclude electrodes at infinity as they dont refer to an entry of the
    % electrode list)
    idx_TX = unique([idx_A; idx_B]);
    idx_TX(idx_TX == 0) = [];
    TX_coo = [all_info.data_block{1}(idx_TX, conf_ele_x), ...
              all_info.data_block{1}(idx_TX, conf_ele_z)];
    n_TX = size(idx_TX , 1);
    
    % Get positions of all electrodes acting as TX.
    idx_RX = unique([idx_M; idx_N]);
    assert(~any(idx_RX == 0), ...
        ['Expect all potential electrodes to be located at known ', ...
        'electrode positions.']);
    RX_coo = [all_info.data_block{1}(idx_RX, conf_ele_x), ...
              all_info.data_block{1}(idx_RX, conf_ele_z)];
    n_RX = size(idx_RX , 1);
          
    % Set TX strengh (i.e. current I).
    % As only rhoa will be considered as supported data type and U is
    % linear in I, the source strengh can be choosen arbitrarily 
    % (e.g. == 1). 
    TX_val = ones(size(TX_coo, 1), 1);
          
    %% Exclude configuration data.
       
    % Check if supported data tokens can be found.
    data_rhoa = strcmp(data_token{2}{5}, all_info.line_header{2}) | ...
                strcmp(data_token_alternative{2}{5}, all_info.line_header{2});
    data_k = strcmp(data_token{2}{9}, all_info.line_header{2}) | ...
                strcmp(data_token_alternative{2}{9}, all_info.line_header{2});
    data_rho = strcmp(data_token{2}{6}, all_info.line_header{2}) | ...
               strcmp(data_token_alternative{2}{6}, all_info.line_header{2});
    data_U = strcmp(data_token{2}{7}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{7}, all_info.line_header{2});
    data_I = strcmp(data_token{2}{8}, all_info.line_header{2}) | ...
             strcmp(data_token_alternative{2}{8}, all_info.line_header{2});
    empty_measure = ~any([data_rhoa, data_rho, data_U, data_I]);
    empty_k = (~any(data_k) || all(all_info.data_block{2}(:,data_k) == 0));
    if any(data_rhoa) && empty_k
        warning(['Data set provides rhoa but no configuration factor k. ', ...'
            'Obtaining k from Neumann formula using actual offsets r ', ...
            'between ekectrodes. If rhoa data is based on a e.g. ', ...
            'topography-corrected k, a systematic error will be introduced.']);
    end

    %% Set up configuration factor k.
    
    if empty_k
        % A           M  N        B
        % |-----------|--|--------|
        %  <--------->             r1
        %  <------------>          r2
        %              <---------> r3
        %                 <------> r4
        % dU = (I rhoa)/2 pi * ((1/r1 - 1/r2) - (1/r3 - 1/r4)) = (I rhoa) / k

        % Obtain RX electrode positions.
        M_coo = [all_info.data_block{1}(idx_M, conf_ele_x), ...
                 all_info.data_block{1}(idx_M, conf_ele_z)];
        N_coo = [all_info.data_block{1}(idx_N, conf_ele_x), ...
                 all_info.data_block{1}(idx_N, conf_ele_z)];

        % Calculate k for respective configutation type.
        getR = @(x) sqrt(x(:,1).^2 + x(:,2).^2);

        % Obtain TX electrodes positions.
        [A_coo, B_coo] = deal(zeros(n_data, 2));
        A_coo(~idx_A_at_inf,:) = ...
            [all_info.data_block{1}(idx_A(~idx_A_at_inf), conf_ele_x), ...
             all_info.data_block{1}(idx_A(~idx_A_at_inf), conf_ele_z)];
        B_coo(~idx_B_at_inf,:) = ...
            [all_info.data_block{1}(idx_B(~idx_B_at_inf), conf_ele_x), ...
             all_info.data_block{1}(idx_B(~idx_B_at_inf), conf_ele_z)];

        % Set position for electrodes at infinity.
        A_coo(idx_A_at_inf,:) = inf(length(find(idx_A_at_inf)), 2);
        B_coo(idx_B_at_inf,:) = inf(length(find(idx_B_at_inf)), 2);

        % Obtain spacings.
        r1 = getR(M_coo - A_coo); % AM
        r2 = getR(N_coo - A_coo); % AN
        r3 = getR(M_coo - B_coo); % BM
        r4 = getR(N_coo - B_coo); % BN

        % Define k.
        k = (2*pi) * 1./((1./r1 - 1./r2) - (1./r3 - 1./r4));
    else
        k = all_info.data_block{2}(:,data_k);
    end
    %% Exclude vector of observed data (=rhoa per definition).
    
    % Obtain rhoa and check uniqueness.
    if ~empty_measure
        % Choose data type to export.
        % (default = rhoa)
        data_var = 'rhoa';
        % For given rhoa information
        if any(data_rhoa)
            rhoa = all_info.data_block{2}(:,data_rhoa);
            % For additonally given rho.
            if any(data_rho) && ~all(all_info.data_block{2}(:,data_rho) == 0)
                if empty_k
                    warning(['Rhoa and rho but no k provided - ', ...
                        'using rho as data']);
                    data_var = 'rho';
                    rho = all_info.data_block{2}(:,data_rho);
                else
                    warning('Rhoa, k and rho provided - ignoring the latter.');
                end
            % For additionally given U and I.
            elseif any(data_U) && any(data_I)
                if ~all(all_info.data_block{2}(:,data_U) == 0) && ...
                        ~all(all_info.data_block{2}(:,data_I) == 0) && ...
                        empty_k
                    warning(['Rhoa and U and I but no k provided - ', ...
                        'using rho as data']);
                    data_var = 'rho';
                    rho = all_info.data_block{2}(:,data_U) ./ ...
                          all_info.data_block{2}(:,data_I);
                elseif ~all(all_info.data_block{2}(:,data_U) == 0) && ...
                           ~all(all_info.data_block{2}(:,data_I) == 0)
                    warning('Rhoa, k and U and I provided - ignoring the two latter.');
                else
                    % Nothing to do.
                end
            end
        % For given rho information
        elseif any(data_rho)
            rhoa = all_info.data_block{2}(:,data_rho) .* k;
            % For additionally given U or I.
            if any(data_U) && any(data_I)
                % TODO: proof consistency.
                warning('Rho and U or I provided - ignoring the latter.');
            end
        % If only U & I are given
        else
            vec_U = all_info.data_block{2}(:,data_U);
            vec_I = all_info.data_block{2}(:,data_I);
            assert(~isempty(length(vec_U)) && (length(vec_U) == length(vec_I)), ...
                'Inconsistency between measured potential and current obtained.');

            rhoa = (all_info.data_block{2}(:,data_U) ./ ...
                    all_info.data_block{2}(:,data_I)) .* k;
        end
    end
    
    %% Set up TX mapping.
    
    % As FE solution is stored as a block vector U = [u1, u2, ...]) 
    % corresponding to the block rhs (i.e. the single TX), handling the 
    % superposition of two point sources can be achieved by superposing the
    % appropriate solution vectors (i.e. two columns).
    % As from the 4-point-formula
    % dU = (I rhoa)/2 pi * ((1/r1 - 1/r2) - (1/r3 - 1/r4))
    % superposition equals to subtract solutions w.r.t. r3 and r4 from 
    % solutions w.r.t. r1 and r2.
    %
    % E.e. A = u2, B = u4 the superposition is given by:
    %     U_FE                  in nDOF x nTX
    % U@allDOF = u2 - u4 = U_FE * [0, 1, 0, -1, 0,...].' = U_FE * r
    % U@allDOF                  in nDOF x nOBS
    %        r =                in  nTX x 1
    %
    % For considering n different combinations of A and B, the nOBS
    % corresponding mapping vectors can be concatenated to a matrix R
    %        R = [r1,r2, ...]   in  nTX x nOBS
    
    % Obtain the indices for electrodes A and B in the list of unique TX
    % positions.
    [~, idx_A_in_idx_TX] = ismember(idx_A, idx_TX);
    [~, idx_B_in_idx_TX] = ismember(idx_B, idx_TX);
    
    % Find combinations where only one pole was used.
    is_A_in_TX_pole_at_B = idx_A_in_idx_TX == 0;
    is_B_in_TX_pole_at_A = idx_B_in_idx_TX == 0;
    is_AB_in_TX_dipole = ~(is_A_in_TX_pole_at_B | is_B_in_TX_pole_at_A);
    idx_AB_in_TX_dipole = find(is_AB_in_TX_dipole);
    n_dipole = length(idx_AB_in_TX_dipole);
    
    % Remove corresponding entries from mappings.
    idx_A_in_idx_TX(idx_A_in_idx_TX == 0) = [];
    idx_B_in_idx_TX(idx_B_in_idx_TX == 0) = [];
        
    % Use sparse matrix definition as R mainly contains zeros.
    % Handle the poles:
    % (As number of entries differ for pole and dipole in i,j, and s,
    % initially these variables are defined as cells 
    % (where for a single entry the number of content does't matter) and
    % finaly are converted to a vector.
    [i, j, s] = deal(cell(n_data, 1));
    % (Here each location is associated with only a single entry in i, j
    % and s).
    s(~is_AB_in_TX_dipole) = {1};
    % Define vector for column indices.
    j(~is_AB_in_TX_dipole) = num2cell(find(~is_AB_in_TX_dipole));
    % Define vector for row indices.
    i(is_A_in_TX_pole_at_B) = num2cell(idx_B_in_idx_TX(is_A_in_TX_pole_at_B));
    i(is_B_in_TX_pole_at_A) = num2cell(idx_A_in_idx_TX(is_B_in_TX_pole_at_A));
    %
    % Handle the dipoles:
    % (Here each location is associated with two entries in i, j
    % and s).
    % s: 1 referring to electrode A, -1 referring to electrode B
    s(is_AB_in_TX_dipole) = mat2cell(kron(ones(n_dipole, 1), [1, -1].'), ...
                                2 + zeros(n_dipole, 1), 1);
    j(is_AB_in_TX_dipole) = mat2cell(kron(idx_AB_in_TX_dipole, [1, 1].'), ...
                                2 + zeros(n_dipole, 1), 1);   
    tmp_row = [idx_A_in_idx_TX(is_AB_in_TX_dipole), ...
               idx_B_in_idx_TX(is_AB_in_TX_dipole)].';
    tmp_row = tmp_row(:);
    i(is_AB_in_TX_dipole) = mat2cell(tmp_row, 2 + zeros(n_dipole, 1), 1);
    
    % Reshape cells to vectors i, j, s.
    s = cell2mat(s);
    j = cell2mat(j);
    i = cell2mat(i);
    
    % Set up sparse matrix.
    R = sparse(i, j, s, n_TX, n_data);
  
    %% Set up RX mapping.
    
    % As the interpolation operator provides the mapping between the FE 
    % solution at every DOF (column index) and the solution at the RX 
    % positions (row index) a subtraction of observations can be obtained
    % by subtracting the corresponding lines of the interpolated FE
    % solution.
    %
    % E.g. M = RX1, N = RX3
    %     U = I * U_FE       in nRX x nTX
    %     I                  in nRX x nDOF
    % dU@RX = [1,0,-1,0,...] * I * U@allDOF
    %       = [1,0,-1,0,...] * I * U * [0,1,0,-1,0,...].'
    %       = l              * I * U * r.'
    
    % Obtain the indices for electrodes M and N in the list of unique RX
    % positions.
    [~, idx_M_in_idx_RX] = ismember(idx_M, idx_RX);
    [~, idx_N_in_idx_RX] = ismember(idx_N, idx_RX);
    
    % Use sparse matrix definition as L mainly contains zeros.
    % Define vector for elements of all l (i.e. L).
    % 1 referring to electrode M, -1 referring to electrode N
    s = kron(ones(n_data, 1), [1, -1].');
    % Define vector for row indices.
    i = kron((1:n_data).', [1, 1].');
    % Define vector for column indices.
    j = [idx_M_in_idx_RX, idx_N_in_idx_RX].';
    j = j(:);
    L = sparse(i, j, s, n_data, n_RX);
    
    %% Summarize.
    
    % Note: to obtain the vector of observed data from the FE block 
    % solution apply: L from left and R from right hand side to the 
    % interpolated solution and exclude the diagonal of the resulting
    % matrix:
    %   dObs = diag(L * (I * U_FE) * R)
    
    dc_conf = struct();
    dc_conf.TX.type = 'point_exact';
    dc_conf.TX.coo = TX_coo;
    dc_conf.TX.val = TX_val;
    dc_conf.RX.coo = RX_coo;
    dc_conf.k = k;
    if ~empty_measure
        switch data_var 
            case {'rhoa'}
                dc_conf.rhoa = rhoa;
            case {'rho'}
                dc_conf.rho = rho;
        end
    end
    dc_conf.mapTX = R;
    dc_conf.mapRX = L;
    if verbosity
        fprintf('done.\n');
    end 
end