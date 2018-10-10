function [fig_handle] = plotRTElement(mesh, fe, varargin)
    % Visualize Raviart-Thomas basis functions on a local/global element.
    %
    % SYNTAX
    %   [fig_handle] = plotRTElement(mesh, fe, idx[, fig_handle])
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %   fe   ... Struct, including all information to set up RT-FE.
    %
    % OPTIONAL PARAMETER
    %   idx  ... Scalar, denoting the cell index for which the basis
    %            functions should be visualized [default = 1].
    %   type ... Character, denoting if cell should be visualized in local
    %            [default] or global cooridnates.
    %   DOF  ... Scalar, global DOF for which corresponding base function
    %            should be plotted.
    %   fig_handle ... Handle to figure, that should be updated.
    %   debugging  ... Logical, use same colors for basis functions
    %                  [default = false]
    %
    % OUTPUT PARAMETER
    %   fig_handle ... Handle to the current figure.
    
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx'})), ...
        'mesh - struct, containing cells info, expected.');
    assert(isstruct(fe) && all(isfield(fe, {'base'})), ...
        'fe - struct, including all information to set up RT-FE, expected.');
    
    % Define possible input keys and its property checks.
    input_keys = {'idx', 'type', 'DOF', 'fig_handle', 'debugging'};
    assertIdx = @(x) assert(isscalar(x), ...
        'idx - scalar, denoting cell index, expected.');
    assertType = @(x) assert(ischar(x) && any(strcmp(x, {'local', 'global'})), ...
        'type - character, denoting coordinate system to be used, expected.');
    assertDOF = @(x) assert(isempty(x) || isscalar(x), ...
        'DOF - scalar, denoting global DOF, expected.');
    assertFig = @(x) assert(isa(x, 'matlab.ui.Figure'), ...
        'fig_handle - figure function handle expected.');
    assertDebug = @(x) assert(islogical(x), ...
        'debugging - logical, true if only one color should be used, expected');
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, 1, assertIdx);
    parser_obj.addParameter(input_keys{2}, 'local', assertType);
    parser_obj.addParameter(input_keys{3}, [], assertDOF);
    parser_obj.addParameter(input_keys{4}, [], assertFig);
    parser_obj.addParameter(input_keys{5}, false, assertDebug);
   
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;   
    
    % Check if desired DOF is related to cell.
    if ~isempty(args.DOF) && ...
            ~ismember(args.DOF, fe.DOF_maps.cell2DOF{args.idx})
         args.DOF = [];
         warning('Provided DOF is not part of current cell.');
    end
    if isempty(args.DOF)
        DOF_plot = 1:3;
    else
        DOF_plot = find(fe.DOF_maps.cell2DOF{args.idx} == args.DOF);
    end
    
    % Create figure if required.
    if isempty(args.fig_handle) || strcmp(args.type, 'local')
        fig_handle = figure();
    else
        fig_handle = args.fig_handle;
    end
    
    %% Create a cell-related grid to plot on.
    
    % local coo.
    tmp1 = num2cell(linspace(0, 1, 5));
    tmp2 = arrayfun(@(x) {ones(1, x)}, 1:length(tmp1));
    y_hat = cellfun(@(x, y) {kron(x, y)}, tmp2, fliplr(tmp1));
    z_hat = arrayfun(@(x) {cell2mat(tmp1(1:x))}, 1:length(tmp1));
    points_hat = [y_hat{:}; z_hat{:}];
    
    % global coo.
    if strcmp(args.type, 'global')
        points = bsxfun(@plus, mesh.maps{args.idx}.B * points_hat, ...
                               mesh.maps{args.idx}.b);
    end
    
    %% Evaluate base functions.
    
    if strcmp(args.type, 'global')
%         trafo = 1/abs(mesh.maps{args.idx}.detB) * mesh.maps{args.idx}.B;
        trafo = 1/(mesh.maps{args.idx}.detB) * mesh.maps{args.idx}.B;
        Phi = arrayfun(@(x, y) {trafo * fe.base.Phi(x, y)}, ...
                  points_hat(1,:), points_hat(2,:)).';
    else
        Phi_hat = arrayfun(@(x, y) {fe.base.Phi(x, y)}, ...
                      points_hat(1,:), points_hat(2,:)).';
    end
    
    %% Plot triangle.
    
    % Add triangle edges;
    hold on
    if strcmp(args.type, 'global')
        patch('Faces', mesh.cell2vtx(args.idx, :), ...
              'Vertices', mesh.vertices, ...
              'FaceColor', 'none', ...
              'EdgeColor', 'blue');
    else
        patch('Faces', [1, 2, 3], ...
              'Vertices', [0, 0; 1, 0; 0, 1], ...
              'FaceColor', 'none');
    end
    hold off
    
    %% Plot base functions.
    
    % Add location of function evaluation.
    hold on
    if strcmp(args.type, 'global')
        plot(points(1,:), points(2,:), '.');
    else
        plot(points_hat(1,:), points_hat(2,:), '.');
    end
    hold off
    
    % Add base function evaluation.
    hold on
    if args.debugging && ~strcmp(args.type, 'local')
        colors = {'black', 'black', 'black'};
    else
        colors = {'yellow', 'black', 'red'}; % colors for base functions
    end
    if strcmp(args.type, 'global')
        % Set length of vectors to be plotted.
        scale = pick(2, 1, 1/7);
        scale = scale * min([max(abs(diff(mesh.cell2cord{args.idx}(:,1)))), ...
                             max(abs(diff(mesh.cell2cord{args.idx}(:,2))))]);
        % Iterate over the locations of function evaluations.
        for ii = 1:size(Phi, 1)
            phi_cur = Phi{ii};
            % Iterate over the three base functions.
            for kk = DOF_plot
                quiver(points(1, ii), points(2, ii), ...
                       phi_cur(1, kk), phi_cur(2, kk), scale, ...
                       'Color', colors{kk}, ...
                       'LineWidth', 1);
            end
        end
    else
        % Set length of vectors to be plotted.
        scale = pick(2, 1, 1/5);
        % Iterate over the locations of function evaluations.
        for ii = 1:size(Phi_hat, 1)
            phi_cur = Phi_hat{ii};
            % Iterate over the three base functions.
            for kk = 1:3
                quiver(points_hat(1, ii), points_hat(2, ii), ...
                       phi_cur(1, kk), phi_cur(2, kk), scale, ...
                       'Color', colors{kk}, ...
                       'LineWidth', 1);
            end
        end
    end
    hold off
end