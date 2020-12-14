function [] = exportObsData2BERT(name, survey, ele, verbosity)
    % Write
    %
    % SYNTAX
    %   [] = exportObsData2BERT(name, survey, ele)
    %
    % INPUT PARAMETER
    %   name   ... Char, denoting the file name to export in.
    %   ele    ... Matrix [n x 2], denoting the electrode positions.
    %   survey ... Struct, containing the tabulated configuration info,
    %              i.e. the combination of A, M, N, B for each observation
    %              and the table column names.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.

    %% Check input.

    assert(ischar(name), ...
        'name - Char, denoting the file name.');
    assert(isstruct(survey) &&  all(isfield(survey, {'data_header', 'data'})), ...
        'survey - Struct, tabulated configuration info and its names, expected.');
    assert(ismatrix(ele) && size(ele, 2) == 2, ...
        'ele - Matrix [n x 2], of electrode positions expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, dennoting if stats should be print, expected.');
    end
    if length(name) < 5 || ~strcmp(name(end-3:end), '.txt')
        name = [name, '.txt'];
    end

    %% Export survey to file in BERT format.

    if verbosity
       fprintf('Export measurement configuration info to BERT file ... ');
    end
    fileID = fopen(name, 'w');
    % Write electrode info block.
    % Set header.
    n_ele = size(ele, 1);
    fprintf(fileID, '%d# Number of electrodes\n', n_ele);
    fprintf(fileID, '# x z\n');
    % Set positions.
    for i = 1:n_ele
        fprintf(fileID, '%.2f %.2f \n', ele(i,:));
    end

    % Write data block.
    % Set header.
    n_data = size(survey.data, 1);
    n_col = size(survey.data, 2);
    fprintf(fileID, '%d# Number of data\n', n_data);
    fprintf(fileID, ['# ', repmat('%s ', 1, n_col), '\n'], survey.data_header{:});
    % Set data.
    for i = 1:n_data
        fprintf(fileID, [repmat('%d ', 1, 4), repmat('%.3f ', 1, n_col-4)], survey.data(i,:));
        if i ~= n_data
            fprintf(fileID, '\n');
        end
    end
    fclose(fileID);
    if verbosity
        fprintf('done.\n');
    end
end
