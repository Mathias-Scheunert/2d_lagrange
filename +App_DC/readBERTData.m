function all_info = readBERTData(name, verbosity)
    % Read in measurement data from BERT formatted file.
    %
    % The unified data format from BERT is defined at
    %   https://gitlab.com/resistivity-net/bert
    %
    % SYNTAX
    %   [dc_conf] = readBERTData(name[, verbosity])
    %
    % INPUT PARAMETER
    %   name ... Char, denoting the file name to read from.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   all_info ... Struct, containing 
    %                - line headers (column names), 
    %                - header units (column name units), 
    %                - the corresponding numerical data for each header
    %                  (stored as matrix of doubles),
    %                - any comments found in the data
    %                To get consistent format for each possible data set,
    %                each field comprises a [3x1] cell array 
    %                 (= number of possible headers in BERT format)
    %                wherein the first dimension of the data_block and
    %                 data_comment as well as the first dimension of the 
    %                 line_header and header_unit arrays always coincide.
    %
    % REMARKS
    %   BERT files are structured as follows:
    %
    %   6# Number of electrodes         <- Header
    %   # x z                           <- Subheader
    %   0     0                         <- Data block line
    %   1     0
    %   ...
    %   6# Number of data
    %   # a b m n rhoa/Ohmm
    %   1   2   3   4  231.2            <- for a pol configuration, the 
    %   ...                                second electrode is set to 0
    %   4# Number of topo points
    %   # x h
    %   0 353.2
    %   ...
    
    %% Check input.
    
    assert(ischar(name), ...
        'name - Char, denoting the BERT formatted file, expected.');
    if nargin < 2
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, dennoting if stats should be print, expected.');
    end

    %% Prepare routine.
    
    % Predefine standardized header lines.
    head_str = {'Number of electrodes', ...
                'Number of data', ...
                'Number of topo points'};
    n_head = length(head_str);
             
    % Predefine standardized data tokens.
    % (These are the 'subheader' column tokens/names).
    data_token = ...
        {{'x', 'y', 'z'}; ...
         {'a',  'b',  'm',  'n',  'rhoa', 'rho', 'err', 'ip', 'ipErr', 'u', 'i', 'k', 'valid'}; ...
         {'x', 'h'} ...
        };
    data_token_alternative = ...
        {{'', '', ''}; ...
         {'c1', 'c2', 'p1', 'p2', 'Ra',   'r',   '',    '',   '',      '',  '',  '',  ''}; ...
         {'', ''} ...
        };
    
    % Predefine meta info tokens.
    comment_tok = '#';
    unit_tok = '/';
    
    %% Read in (text)file and separate.
    
    if verbosity
        fprintf('Read out data from BERT formatted file ... ');
    end
    % Get whole content from file.
    file_content = fileread(name);
    
    % Replace ASCII formatted tabs.
    % TODO: check if there might be further strange formatting issues to
    %       handle here.
    file_content = regexprep(file_content, char(9), ' ');

    
    % Store each line of the input file in a seperate cell.    
    file_content = textscan(file_content,  '%s', 'Delimiter', '\n');
    file_content = file_content{1};
    
    % Delet leading and tailing whitspaces.
    file_content = strtrim(file_content);
    
    % Check file.
    assert(~isempty(file_content), 'Empty or corrupt file detected.');
    n_line = size(file_content, 1);
        
    %% Go through content and search for line headers.
    
    head_line_idx = false(n_line, n_head);
    for ii = 1:n_head
        % Find header line index by searching for token.
        head_line_idx(:,ii) = contains(file_content, head_str{ii}, ...
                                       'IgnoreCase', true);
                                   
        % If no token was provided, search for the subheader lines,
        % starting with a comment, instead.
        if all(~head_line_idx(:,ii))           
           % Search for comment token in every line.
           if ~exist('subhead_line_idx', 'var')
               subhead_line_idx = cellfun(@(x) strcmp(x(1), comment_tok), file_content);
               subhead_line_idx = find(subhead_line_idx);
               n_subhead_line = length(subhead_line_idx);
               
               % Check consistency to expected BERT format.
               assert(n_subhead_line == 2 || n_subhead_line == 3, ...
                   ['Expect file to contain at least two information ', ...
                   'blocks. That are positions of electrodes and the ', ...
                   'configuration/measurement data.']);
           end
           
           % Set current headline index.
           % (Only where subheader line indices where found)
           if ii <= n_subhead_line
            head_line_idx(subhead_line_idx(ii) - 1, ii) = true;
           end
        end
    end
    
    %% Read out header information.

    head_info = cell(n_head, 1);
    [subhead_info, subhead_unit_info] = deal(cell(n_head, 1));
    for ii = 1:n_head
        cur_line_idx = find(head_line_idx(:,ii));
        if isempty(cur_line_idx)
            % Skip case as this header doesn't exist.
        else
            % Exclude numeric header content (i.e. number of respective 
            % entries) by reading string up to comment token.
            head_info{ii} = cell2mat(...
                textscan(file_content{cur_line_idx}, ...
                         ['%d[^', comment_tok,']']));

            % Exclude column names from subheader (including units).
            cur_subhead = regexprep(file_content{cur_line_idx + 1}, ...
                                    comment_tok, '');
            subhead_info{ii} = textscan(cur_subhead, '%s');
            subhead_info{ii} = subhead_info{ii}{:}.';
            
            % Exclude units from column names.
            subhead_unit_info{ii} = cell(size(subhead_info{ii}));
            unit_defined = contains(subhead_info{ii}, unit_tok);
            if any(unit_defined)
                % Separate names and units.
                tmp_units = cellfun(@(x) {strsplit(x, unit_tok)}, ...
                    subhead_info{ii}(unit_defined)).';
                tmp_units = vertcat(tmp_units{:});
                
                % Just keep names.
                subhead_info{ii}(unit_defined) = tmp_units(:,1);
                
                % Store corresponding units.
                subhead_unit_info{ii}(unit_defined) = tmp_units(:,2);
            end
            
            % Check for consistency (case insensitive).
            % I.e. check if all obtained column names are part of the 
            % predefined/known (alternative) data tokens.
            cur_token = data_token{ii};
            cur_alt_token = data_token_alternative{ii};
            token_in_subhead = cellfun(@(x) ...
                {strcmpi(cur_token, x)}, subhead_info{ii});
            alt_token_in_subhead = cellfun(@(x) ...
                {strcmpi(cur_alt_token, x)}, subhead_info{ii});
            % TODO: may weaken this assert and just give a warning here.
            assert(all(cellfun(@(x, y) any(x) || any(y), ...
                token_in_subhead, alt_token_in_subhead)), ...
                sprintf(['Some column names of "%s" info are not part of ', ...
                'list of known names.'], head_str{ii}));
        end
    end
    
    %% Read out corresponding data for each header block.
    
    % Get info blocks as cell array of strings.
    % TODO: add check, proofing if number of data lines given in header
    %       really fits to the number of lines following this header.
    %       E.g. use the head_line_idx of the next header and compare it
    %       with the current head_line_idx (+1) + head_info.
    [block_info, block_comment] = deal(cell(n_head, 1));
    for ii = 1:n_head
        cur_start_line = find(head_line_idx(:,ii)) + 2;
        cur_length = head_info{ii} - 1;
        block_info{ii} = file_content(cur_start_line:(cur_start_line + cur_length));
        block_comment{ii} = cell(size(block_info{ii}));
    end
    
    % Exclude data from stings in each cell.
    empty_block = cellfun(@isempty, block_info);
    for ii = 1:n_head
        if empty_block(ii)
            % Skip empty block.
        else
            % Allocate data matrix from number of lines and column names.
            n_cur_line = size(block_info{ii}, 1);
            n_cur_col = size(subhead_info{ii}, 2);
            cur_data_mat = zeros(n_cur_line, n_cur_col);
            
            % Exclude numerical data and comment string.
            % Loop through lines.
            % TODO: may speed up.
            %       -> As we expect only a moderate number of lines, this 
            %       brute force approach might be ok for now.
            for ll = 1:n_cur_line
                % Separate the entries of the current line.
                cur_line_split = strsplit(block_info{ii}{ll});
                
                % Store the comment.
                line_commment = contains(cur_line_split, comment_tok);
                if any(line_commment)
                    idx_comment_start = find(line_commment);
                    assert(idx_comment_start > n_cur_col, ...
                        ['Comments expected at the end of the data. ', ...
                        'Please check file format.']);
                    % Merge text cells to string (exclude comment token).
                    tmp_comment = cur_line_split(idx_comment_start:end);
                    tmp_comment(1:end-1) = cellfun(@(x) ...
                                        {[x, ' ']} ,tmp_comment(1:end-1));
                    block_comment{ii}{ll} = [tmp_comment{:}];
                    block_comment{ii}{ll} = block_comment{ii}{ll}(2:end);
                    
                    % Remove comment cells.
                    cur_line_split = cur_line_split(1:n_cur_col);
                end
                
                % Store the numerical data.
                cur_data_mat(ll,:) = cellfun(@str2num, cur_line_split);
            end
            
            % Replace block info.
            block_info{ii} = cur_data_mat; 
        end
    end
    
    %% Check consistency of topography info.
    
    % Check if relevant data was obtained.
    warn_text = cell(0);
    if isempty(block_info{1}) || isempty(block_info{2})
        cur_text = 'Electrode or configuration data was missing in file.';
        warning(cur_text);
        warn_text = [warn_text; cur_text];
    else 
        dim = num2str(length(subhead_info{1}));
        type = [dim, 'D'];
    end
    
    if ~isempty(block_info{3}) && ~isempty(block_info{1})
       % Check if topography (x-)location matches electrode positions.
       x_in_ele = strcmp(subhead_info{1}, 'x');
       z_in_ele = strcmp(subhead_info{1}, 'z');
       topo_x = block_info{3}(:,1);
       topo_z = block_info{3}(:,2);
       [~, ele_in_topo] = ismember(topo_x, block_info{1}(:,x_in_ele));
       if any(ele_in_topo == 0)
           cur_text = ['Some topography locations are not part of ', ...
               'electrode location list.'];
           warning(cur_text);
           warn_text = [warn_text; cur_text];
           ele_in_topo(ele_in_topo == 0) = [];
       end
       
       % Check if topography (z-)location matches electrode positions.
       ele_z = block_info{1}(ele_in_topo, z_in_ele);
       match_z = topo_z(ele_z ~= 0) == ele_z(ele_z ~= 0);
       if ~all(match_z)
           cur_text = ['Some heights of electrodes do not match the ', ...
               'topography info at this location.'];
           warning(cur_text);
           warn_text = [warn_text; cur_text];
       end
    end
    
    %% Summarize output.
    
    all_info = struct();
    if exist('type', 'var')
        all_info.type = type;
    end
    all_info.line_header = subhead_info;
    all_info.header_unit = subhead_unit_info;
    all_info.data_block = block_info;
    all_info.data_comment = block_comment;
    if ~isempty(warn_text)
        all_info.warnings = warn_text;
    end
    if verbosity
        fprintf('done.\n');
    end 
end