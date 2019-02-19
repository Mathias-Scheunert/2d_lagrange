function [] = plotGradient(fe, mesh, u, varargin)
    % Visualize the normalized gradient of the solution at cell centroids.
    %
    % The gradients are evaluated at cells centroids but plotted at cells
    % incenter!
    %
    % Gradient of the solution at point of the cell k is given by:
    %   \grad Phi_k = \sum_i(k) ( B_k^(-T) \grad(\phi_i(k)(x_mid)) * u_i(k))
    % for
    %   k = index of cell
    %   i = index of DOF of cell k
    %
    % SYNTAX
    %   [] = plotGradient(fe, mesh, u[, varargin])
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %   u    ... Vector, providing the solution of the FE fwd problem at
    %            the DOF.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed. [default = false]
    %   param     ... Vector [n x 1], denoting the cell parameters.
    %                 [default = NaN*[number_of_cells x 1]]
    %   sign      ... Character, denoting if gradient is oriented on 'pos'
    %                 or 'neg' direction. [default = 'pos']
    %
    % REMARKS
    %   Note that the length of the plotted gradient vector arrows in a 
    %   cell does not represent the amplitude of the gradient. It is rather
    %   derived from the cell's inner cirle diameter!

    %% Check input
    
    assert(isstruct(fe) && all(isfield(fe, {'DOF_maps', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2cord', 'bnd_edge'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isvector(u) && length(u) == fe.sizes.DOF, ...
        'u - vector, containing the solution for at all DOF, expected.');
    
    % Define optional input keys and its properties checks.
    input_keys = {'param', 'verbosity', 'sign'};
    assertParam = @(x) assert(isvector(x) && length(x) == fe.sizes.cell, ...
        'param - vector, including the constant parameter values, expected');
    assertVerbose = @(x) assert(islogical(x), ...
        'verbosity - logical, denoting if status should be printed, expected');
    assertSign = @(x) assert(ischar(x) && any(strcmp(x, {'neg', 'pos'})), ...
        'sign - character, denoting if pos. or neg grad is plotted, expected');
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, NaN*zeros(fe.sizes.cell, 1), assertParam);
    parser_obj.addParameter(input_keys{2}, false, assertVerbose);
    parser_obj.addParameter(input_keys{3}, 'pos', assertSign);
    
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    % Fetch info.
    n_cell = fe.sizes.cell;
    
    %% Evaluate solution at cell centroids.
    if isnan(args.param)
        % Get local basis functions at cell centroids.
        loc_basis = fe.base.Phi(1/3, 1/3);

        % Evaluate solution at all cell centroids.
        u_cell = zeros(n_cell, 1);
        for ii = 1:n_cell                
            % Get current DOF.
            cur_DOF = fe.DOF_maps.cell2DOF{ii};

            % Add DOF and sum up solution.
            u_cell(ii) = sum(bsxfun(@times, u(cur_DOF).', loc_basis), 2);
        end
    end
        
    %% Get gradient of solution.
    
    if args.verbosity
       fprintf('Evaluate gradient of the solution ... '); 
    end
    % Evaluate gradients of basis function at these points and superpose 
    % evaluations for all basis functions of a cell.
    grad_coo = zeros(n_cell, 2);
    cell_incenter = zeros(n_cell, 2);
    
    % Evaluate gradient of basis functions at cell centroid in
    % barycentric coordinates.
    loc_base = fe.base.grad_Phi(1/3, 1/3);
    
    % Loop over all cells.
    for ii = 1:n_cell
        % Get DOF for current cell.
        cur_DOF = fe.DOF_maps.cell2DOF{ii};
                
        % Map to global coordinates.
        cur_glo_base = mesh.maps{ii}.BinvT * loc_base;
        
        % Add DOF and sum up.
        cur_glo_grad =  sum(bsxfun(@times, u(cur_DOF).', cur_glo_base), 2);
        
        % Normalize.
        grad_coo(ii,:) = cur_glo_grad / norm(cur_glo_grad);
        
        % Get current cell's incircle (see wikipedia).
        % Get edge length's.
        cur_edges = mesh.cell2edg(ii,:);
        cur_edges_cords = mesh.edge2cord(cur_edges);
        cur_edges_length = cell2mat(cellfun(@(x) ...
            {sqrt((x(2,1) - x(1,1))^2 + (x(2,2) - x(1,2))^2)}, ...
            cur_edges_cords));
        % Get semiperimeter.
        cur_s = sum(cur_edges_length)/2;
        % Get inner circle radius.
        cur_r = sqrt(prod(cur_s - cur_edges_length)/cur_s);
        
        % Get cell incenter in global coordinates.
        cell_incenter(ii, :) = sum(flipud(cur_edges_length) .* ...
                                   mesh.cell2cord{ii}, 1) / ...
                                      (cur_s * 2);
        
        % Adjust length of the gradient to be little smaler than the 
        % current cells inner circle diameter.
        grad_coo(ii,:) = grad_coo(ii,:) * 1.8*cur_r;
        switch args.sign
            case 'pos'
                % Nothing to do.
            case 'neg'
                grad_coo(ii,:) = -grad_coo(ii,:);
        end
    end
    if args.verbosity
       fprintf('done.\n'); 
    end
    
    %% Visualize gradient.
    
    if args.verbosity
       fprintf('Print gradient of the solution ... '); 
    end
    % Create figure.
    figure();
    axis('equal');
    set(gca, 'Ydir', 'reverse') % As y should point downwards.
    xlim([min(mesh.vertices(:,1)), max(mesh.vertices(:,1))]);
    ylim([min(mesh.vertices(:,2)), max(mesh.vertices(:,2))]);
    caxis('auto');

    % Visualize solution with underlying mesh.
    hold on
    if isnan(args.param)
        patch('Faces', mesh.cell2vtx, 'Vertices', mesh.vertices, ...
              'FaceVertexCData', u_cell, ...
              'FaceColor', 'flat');
        h = colorbar();
        ylabel(h, 'potential');
    else
        patch('Faces', mesh.cell2vtx, 'Vertices', mesh.vertices, ...
              'FaceVertexCData', args.param, ...
              'FaceColor', 'flat');
        h = colorbar();
        ylabel(h, 'cell parameter value');
    end
    hold off
    
    % Slightly shift them, sucht that center of plotted vectors will be
    % located at the cell incenter.
    % (otherwise the vector origin will be there)
    plot_mid = [cell_incenter(:,1) - grad_coo(:,1)./2, ...
                cell_incenter(:,2) - grad_coo(:,2)./2];
    
    % Display the solution for all DOF (at an adapted triangulation).
    hold on
    quiver(plot_mid(:,1), plot_mid(:,2), grad_coo(:,1), grad_coo(:,2), ...
        'Color', 'red', ...
        'LineWidth', 2, ...
        'AutoScale','off');
    hold off
    
    % Adjust figure.
    xlabel('x');
    ylabel('y');
    if args.verbosity
       fprintf('done.\n'); 
    end
end