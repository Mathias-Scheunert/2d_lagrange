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
    %   cell_idx   ... Scalar, denoting the cell index for which the basis
    %                  functions should be visualized [default = 1].
    %   type       ... Character, denoting if cell should be visualized in
    %                  local [default] or global coordinates.
    %   DOF        ... Scalar, global DOF for which corresponding base
    %                  function should be plotted.
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
    input_keys = {'cell_idx', 'type', 'DOF', 'fig_handle', 'debugging'};
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
            ~ismember(args.DOF, fe.DOF_maps.cell2DOF{args.cell_idx})
         args.DOF = [];
         warning('Provided DOF is not part of current cell.');
    end
    if isempty(args.DOF)
        DOF_loc_plot = 1:3;
        DOF_glo_plot = mesh.loc2glo(DOF_loc_plot);
    else
        % Note: As cell2DOF contains the loc2glo-mapping, find() will
        % provide the local DOF ordering (see FeRT.getDOFMap).
        DOF_loc_plot = find(fe.DOF_maps.cell2DOF{args.cell_idx} == args.DOF);
        % I.e. here the mesh.cell2edge ordering will be obtained.
        DOF_glo_plot = mesh.loc2glo(DOF_loc_plot);
    end

    % Create figure if required.
    if isempty(args.fig_handle) || strcmp(args.type, 'local')
        fig_handle = figure();
    else
        fig_handle = args.fig_handle;
    end

    %% Create a cell-related grid to plot on.

    % local coo.
    tmp1 = linspace(0, 1, 5);
    if ~isscalar(args.DOF)
        % Distributed over the triangle.
        tmp1 = num2cell(tmp1);
        tmp2 = arrayfun(@(x) {ones(1, x)}, 1:length(tmp1));
        y_hat = cellfun(@(x, y) {kron(x, y)}, tmp2, fliplr(tmp1));
        z_hat = arrayfun(@(x) {cell2mat(tmp1(1:x))}, 1:length(tmp1));
        points_hat = [y_hat{:}; z_hat{:}];
    else
        % Only at the considered edge.
        tmp2 = {[tmp1; 0*tmp1], [tmp1; fliplr(tmp1)], [0*tmp1; tmp1]};
        points_hat = tmp2{DOF_loc_plot};
    end

    % global coo.
    if strcmp(args.type, 'global')
        points = bsxfun(@plus, mesh.maps{args.cell_idx}.B * points_hat, ...
                               mesh.maps{args.cell_idx}.b);
    end

    %% Evaluate base functions.

    if strcmp(args.type, 'global')
        % Get all edges global normal of current cell (in global ordering).
        glo_edge_normals = Mesh.getEdgeNormal(mesh, 1:size(mesh.edge2vtx, 1));
        glo_edge_normals = cell2mat(cellfun(@(x) {x{1}}, glo_edge_normals));
        cur_cell_global_normals = glo_edge_normals(...
                                      mesh.cell2edg(args.cell_idx,:), :);
        % Get current edge local normal (in global ordering).
        cur_cell_normals = Mesh.getEdgeNormal(mesh, ...
                               mesh.cell2edg(args.cell_idx,:), args.cell_idx);
        cur_cell_normals = cell2mat([cur_cell_normals{:}].');

        % Get edge indices to search for.
        cur_edge_idx = mesh.cell2edg(args.cell_idx,:);

        % Search for occurence of current edge(s) in cell2edge list.
        edge_1_occur = mesh.cell2edg == cur_edge_idx(1);
        edge_2_occur = mesh.cell2edg == cur_edge_idx(2);
        edge_3_occur = mesh.cell2edg == cur_edge_idx(3);

        % Obtain indes for first occurence.
        [edge_1_2_cell, ~] = find(edge_1_occur, 1, 'first');
        [edge_2_2_cell, ~] = find(edge_2_occur, 1, 'first');
        [edge_3_2_cell, ~] = find(edge_3_occur, 1, 'first');

        % Derive sign from position w.r.t current cell index.
        Phi_sign = [2 * (edge_1_2_cell == args.cell_idx) - 1, ...
                    2 * (edge_2_2_cell == args.cell_idx) - 1, ...
                    2 * (edge_3_2_cell == args.cell_idx) - 1];
        % -> as Phi_sign will act on the local basis functions no mapping
        % is incorporated here (see FeRT.assembleMassDiv)

        % By using abs(), the basis functions are always pointing outwards.
        % (Note the local ordering here!).
        trafo = 1/abs(mesh.maps{args.cell_idx}.detB) * mesh.maps{args.cell_idx}.B;
        Phi = arrayfun(@(x, y) {trafo * fe.base.Phi(x, y)}, ...
                  points_hat(1,:), points_hat(2,:)).';

        % Correct orientation and ordering.
        % (Note reordering from local to global via mesh.glo2loc!)
        Phi = cellfun(@(x) {x(:,mesh.glo2loc) .* Phi_sign}, Phi);

        % Proof normal component.
        fprintf(sprintf(['Cell: %d Phi: ', ...
                        repmat('%d, ', 1, length(DOF_loc_plot)), ...
                        '(local ordering)\n'], ...
                args.cell_idx, DOF_loc_plot));
        if isscalar(args.DOF)
            % Get normal vector of current edge in global coords.
            cur_edge_normal = Mesh.getEdgeNormal(mesh, args.DOF, args.cell_idx);
            cur_edge_normal = cur_edge_normal{1}{1};
            % Check:
            assert(isequal(cur_cell_normals(DOF_glo_plot,:), ...
                           cur_edge_normal), ...
                   'Inconsistent normals - check ordering of DOF.');

            % Get basis function evaluation at edge DOF in global coords.
            % (Note the local ordering here!).
            Phi_at_DOF = trafo * ...
                fe.base.Phi(fe.base.DOF(DOF_loc_plot,1), ...
                            fe.base.DOF(DOF_loc_plot,2));
            % Select apprpriate Phi.
            Phi_at_DOF = Phi_at_DOF(:, mesh.glo2loc);
            Phi_at_DOF = Phi_at_DOF(:, DOF_glo_plot);

            % Check if basis function points in the direction of the
            % current edge normal.
            % (should always be the case, if abs() is used in trafo!)
            fprintf(sprintf('Edge normal points in direction of Phi: %d\n', ...
                dot(Phi_at_DOF.', cur_edge_normal) > 0));

            % Check if basis function points in the direction of the global
            % normal.
            % (should only be the case for one of the two adjacent cells!)
            fprintf(sprintf(['Edge normal w.r.t. cell %d points in ', ...
                'direction of the global edge normal: %d\n'], ...
                args.cell_idx, ...
                dot(cur_cell_global_normals(DOF_glo_plot,:), cur_edge_normal) > 0));

%             % is equal to (if abs() is used)
%             dot(cur_cell_global_normals(DOF_glo_plot,:), Phi_at_DOF.') > 0
        end
    else
        Phi_hat = arrayfun(@(x, y) {fe.base.Phi(x, y)}, ...
                      points_hat(1,:), points_hat(2,:)).';
    end

    %% Process legend, if available.

    leg_handl = findobj(fig_handle, 'Type', 'Legend');
    if ~isempty(leg_handl)
        leg_handl_str = leg_handl.String;
        leg_handl.AutoUpdate = 'off';
    end

    %% Plot triangle.

    % Add triangle edges;
    hold on
    if strcmp(args.type, 'global')
        patch('Faces', mesh.cell2vtx(args.cell_idx, :), ...
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
    colors = {'yellow', 'black', 'red'}; % colors for base functions
    if strcmp(args.type, 'global')
        % Set length of vectors to be plotted.
        scale = pick(2, 1, 1/7);
        scale = scale * min([max(abs(diff(mesh.cell2cord{args.cell_idx}(:,1)))), ...
                             max(abs(diff(mesh.cell2cord{args.cell_idx}(:,2))))]);
        % Iterate over the locations of function evaluations.
        plot_handl = cell(length(DOF_glo_plot), 1);
        plot_names = {'Phi_{1, glob}', 'Phi_{2, glob}', 'Phi_{3, glob}'};
        for ii = 1:size(Phi, 1)
            phi_cur = Phi{ii};
            % Iterate over the three base functions.
            if ~isempty(leg_handl) && ii == size(Phi, 1)
                % Allow to update legend
                leg_handl.AutoUpdate = 'on';
            end
            for kk = DOF_glo_plot
                plot_handl{kk} = quiver(points(1, ii), points(2, ii), ...
                       phi_cur(1, kk), phi_cur(2, kk), scale, ...
                       'Color', colors{kk}, ...
                       'LineWidth', 2);
            end
            if kk == DOF_glo_plot(end) && ii == size(Phi, 1)
                leg_handl = findobj(fig_handle, 'Type', 'Legend');
                if isempty(leg_handl)
                    % Create new legend.
                    legend([plot_handl{DOF_glo_plot}], ...
                            plot_names{DOF_glo_plot});
                else
                    % Add new strings to the updated legend.
                    leg_handl.String = [leg_handl_str, plot_names{DOF_glo_plot}];
                end

            end
        end
    else
        % Set length of vectors to be plotted.
        scale = pick(2, 1, 1/5);
        % Iterate over the locations of function evaluations.
        plot_handl = cell(3, 1);
        for ii = 1:size(Phi_hat, 1)
            phi_cur = Phi_hat{ii};
            % Iterate over the three base functions.
            for kk = 1:3
                plot_handl{kk} = quiver(points_hat(1, ii), points_hat(2, ii), ...
                       phi_cur(1, kk), phi_cur(2, kk), scale, ...
                       'Color', colors{kk}, ...
                       'LineWidth', 2);
                if kk == 3 && ii == size(Phi_hat, 1)
                    legend([plot_handl{:}], ...
                        'Phi_{1, loc}', 'Phi_{2, loc}', 'Phi_{3, loc}');
                end
            end
        end
    end
    hold off
end
