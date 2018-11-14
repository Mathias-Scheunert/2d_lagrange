function [] = plotSolution(fe, mesh, u, varargin)
    % Visualize the solution of the 2D FE fwd problem.
    % 
    % SYNTAX
    %   [] = plotSolution(fe, mesh, u[, varargin])
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
    %   param     ... Vector, including the constant parameter value for
    %                 each cell. [default = false]
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.   [default = false]
    %   debug     ... Logical, denoting if additional info to debug will be
    %                 printed.   [default = false]
    %   style     ... Character, denoting if solution will be plotted in
    %                 '2D' or '3D'. [default = '2D']
    
    %% Check input
    
    assert(isstruct(fe) && all(isfield(fe, {'DOF_maps', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2cord', 'bnd_edge'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isvector(u) && length(u) == fe.sizes.DOF, ...
        'u - vector, containing the solution for at all DOF, expected.');
    
    % Define optional input keys and its properties checks.
    input_keys = {'param', 'verbosity', 'debug', 'style'};
    assertParam = @(x) assert(isvector(x) && length(x) == fe.sizes.cell, ...
        'param - vector, including the constant parameter values, expected');
    assertVerbose = @(x) assert(islogical(x), ...
        'verbosity - logical, denoting if status should be printed, expected');
    assertDebug = @(x) assert(islogical(x), ...
        'debug - logical, denoting if debug stuff is printed, expected');
    assertStyle = @(x) assert(ischar(x) && any(strcmp(x, {'2D', '3D'})), ...
        'style - character, denoting if figure is plotted in 2D or 3D, expected');
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, zeros(fe.sizes.cell, 1), assertParam);
    parser_obj.addParameter(input_keys{2}, false, assertVerbose);
    parser_obj.addParameter(input_keys{3}, false, assertDebug);
    parser_obj.addParameter(input_keys{4}, '2D', assertStyle);
   
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;
    
    %% Get coordinates for all DOF.

    x = fe.DOF_maps.DOF_coo(:,1);
    y = fe.DOF_maps.DOF_coo(:,2);
    if args.debug
        
        if args.verbosity
           fprintf('Print mesh and nodes ... '); 
        end
        % Plot mesh.
        Plot.plotMesh(mesh);
        hold on
          % Add DOF locations.
          scatter(x(1:fe.sizes.vtx), y(1:fe.sizes.vtx), 14, 'blue');
          if fe.order == 2
              scatter(x(fe.sizes.vtx + 1:end), y(fe.sizes.vtx + 1:end), 26, 'red');
          else
              error('Visualization of order higher than two not supported yet.');
          end
        hold off
        
        % Plot mesh.
        hold on
        cell_coo_all = cell2mat(mesh.cell2cord);
        cell_coo_x = reshape(cell_coo_all(:,1), [3, fe.sizes.cell]);
        cell_coo_y = reshape(cell_coo_all(:,2), [3, fe.sizes.cell]);
        patch(cell_coo_x, cell_coo_y, max(u) * (param / max(param)));
        hold off
        if args.verbosity
           fprintf('done.\n'); 
        end
    end
    
    %% Add solution.
    
    if args.verbosity
       fprintf('Print solution ... '); 
    end
    figure();
    axis('equal');
    set(gca, 'Ydir', 'reverse') % As y should point downwards.
    xlim([min(mesh.vertices(:,1)), max(mesh.vertices(:,1))]);
    ylim([min(mesh.vertices(:,2)), max(mesh.vertices(:,2))]);
    caxis('auto');
    % Display the solution for all DOF (at an adapted triangulation).
    switch fe.order
        case 1
            if strcmp(args.style, '3D')
                trisurf(mesh.cell2vtx, x, y, u, 'edgecolor', 'none');
            else
                patch('Faces', mesh.cell2vtx, 'Vertices', mesh.vertices, ...
                      'FaceVertexCData', u, ...
                      'FaceColor', 'flat', 'edgecolor', 'none');
            end
        case 2
            mesh_plot = mesh;
            mesh_plot.type = 'basic';
            mesh_plot = Mesh.refineMeshUniform(mesh_plot, 1);
            [~, DOF2vtx_map] = ismember(mesh_plot.vertices, ...
                fe.DOF_maps.DOF_coo, 'rows');
            if strcmp(args.style, '3D')
                trisurf(mesh_plot.cell2vtx, ...
                    x(DOF2vtx_map), y(DOF2vtx_map), u(DOF2vtx_map), ...
                    'edgecolor', 'none');
            else
                patch('Faces', mesh_plot.cell2vtx, ...
                      'Vertices', mesh_plot.vertices, ...
                      'FaceVertexCData', u(DOF2vtx_map), ...
                      'FaceColor', 'flat', 'edgecolor', 'none');
            end
    end
    xlabel('x');
    ylabel('y');
    h = colorbar();
    ylabel(h, 'solution amplitude');
    if args.verbosity
       fprintf('done.\n'); 
    end
end