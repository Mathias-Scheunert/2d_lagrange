function dc_conf = createConfigBERT(ele, type, varargin)
    % Create DC measurement konfiguration data in BERT format.
    %
    % SYNTAX
    %   dc_conf = createConfigBERT(ele, type, varargin)
    %
    % INPUT PARAMETER
    %   ele  ... Matrix [n x 2], denoting the electrode positions.
    %   type ... Char, denoting the desired configuration type.
    %
    % OPTIONAL PARAMETER
    %   store     ... Logical, denoting if configuration will be stored in
    %                 BERT format. [default = false]
    %   name      ... Char, denoting file name to be stored in. 
    %                 [default = 'configuration']
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
    %               - configuration factor k.
    
    %% Check input.
    
    assert(ismatrix(ele) && size(ele, 2) == 2, ...
        'ele - Matrix [n x 2], denoting the electrode positions, expected.');
    assert(ischar(type) && any(strcmp(type, {'wenner', 'schlumberger', ...
                                    'pol-dipol-fw', 'pol-dipol-rw'})), ...
        ['type - Char, denoting the DC measurement configuration, ', ...
         'expected.']);
    
    % Define possible input keys and its properties checks.
    input_keys = {'store', 'name', 'verbosity'};
    assertStore = @(x) assert(islogical(x), ...
        'store - Logical, denoting if configuration should be stored in .dat file');
    assertName = @(x) assert(ischar(x), ...
        'name - Char for file name to be stored, expected.');
    assertVerbose = @(x) assert(islogical(x), ...
        'verbosity - logical, denoting if status should be printed, expected');
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, false, assertStore);
    parser_obj.addParameter(input_keys{2}, 'configuration.txt', assertName);
    parser_obj.addParameter(input_keys{3}, false, assertVerbose);
   
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    %% Construct survey.

    if args.verbosity
       fprintf('Construct measurement configuration info ... '); 
    end 
    n_ele = size(ele, 1);
    switch type
        case 'wenner'
            % General scheme:
            %   A-ka-M-ka-N-ka-B
            % Moving along the electrode chain for variing n, e.g.:
            % |-a-|-a-|-a-|-a-|-a-|-a-|-a-|-a-|-a-|
            % 1a:
            % A-a-M-a-N-a-B---|---|---|---|---|---|
            % |---A-a-M-a-N-a-B---|---|---|---|---|
            % ...
            % 2a:
            % A--2a---M--2a---N--2a---B---|---|---|
            % |---A--2a---M--2a---N--2a---B---|---|
            % 3a:
            % A-----3a----M-----3a----N-----3a----B
            
            % Get max spacing factor k:
            k_max = floor((n_ele - 1) / 3);
            
            % Get max number of measurements.
            n_active_ele = 4;
            n_obs_max = sum((n_ele - (n_active_ele - 1) * (1:k_max)));
            
            % Define matrix of electrode combinations.
            survey = struct();
            survey.data = zeros(n_obs_max, n_active_ele);
            survey.data_header = {'A', 'M', 'N', 'B'};
            cur_row_idx = 0;
            for ii = 1:k_max
                % Get number of observations for current spacing.
                cur_n_obs = n_ele - (n_active_ele - 1) * ii;
                
                % Fill matrix rows.
                for jj = 1:cur_n_obs
                    % Row entries are: [jj, jj + ii,  jj + 2*ii,  jj + 3*ii]
                    survey.data(cur_row_idx + jj,:) = ...
                        jj + (((1:n_active_ele)  - 1) * ii);
                end
                
                % Update auxiliary variable.
                cur_row_idx = cur_row_idx + cur_n_obs;
            end
        otherwise 
            error('Not implemented yet');
    end
    
    %% (Temporarily) export information.
    
    % Make sure that the .dat ending is included in name string.
    if ~strcmp(args.name(end - 3:end), '.txt')
        args.name = [args.name, '.txt'];
    end
    App_DC.exportObsData2BERT(args.name, survey, ele);
    
    %% Read in complete configuration information.

    info = App_DC.getConfigFromBERTData(args.name);
    if ~args.store
        delete(args.name);
    end
    
    %% Summarize.
    
    dc_conf = struct();
    dc_conf.TX = info.TX;
    dc_conf.RX = info.RX;
    dc_conf.k = info.k;
    dc_conf.mapTX = info.mapTX;
    dc_conf.mapRX = info.mapRX;
    if args.verbosity
        fprintf('done.\n');
    end 
end