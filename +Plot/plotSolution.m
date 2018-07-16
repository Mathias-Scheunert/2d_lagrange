function [] = plotSolution(fe, mesh, u, param, verbosity)
    % Visualize the solution of the 2D FE fwd problem.
    % 
    % SYNTAX
    %   [] = plotSolution(fe, mesh, u, param, verbosity)
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %   u    ... Vector, providing the solution of the FE fwd problem at
    %            the DOF.
    %
    % OPTIONAL PARAMETER
    %   param     ... Vector, including the constant parameter value for
    %                 each cell.
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    
    debug = pick(1, false, true);
    
    %% Check input
    
    assert(isstruct(fe) && all(isfield(fe, {'DOF_maps', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2cord', 'bnd_edge'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isvector(u) && length(u) == fe.sizes.DOF, ...
        'u - vector, containing the solution for at all DOF, expected.');
    if nargin < 5
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    if nargin < 4
        param = zeros(fe.sizes.cell, 1);
    else
        assert(isvector(param) && length(param) == fe.sizes.cell, ...
            'param - vector, including the constant parameter values, expected');
    end
    
    %% Get coordinates for all DOF.

    x = fe.DOF_maps.DOF_coo(:,1);
    y = fe.DOF_maps.DOF_coo(:,2);
    if debug
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
        pause();
        close(gcf);
    end
    
    %% Create figure.
    
    % Show the underlying grid.
    figure();
    set(gca, 'Ydir', 'reverse')
    xlim([min(mesh.vertices(:,1)), max(mesh.vertices(:,1))]);
    ylim([min(mesh.vertices(:,2)), max(mesh.vertices(:,2))]);

    if verbosity
       fprintf('Print mesh ... '); 
    end
    cell_coo_all = cell2mat(mesh.cell2cord);
    cell_coo_x = reshape(cell_coo_all(:,1), [3, fe.sizes.cell]);
    cell_coo_y = reshape(cell_coo_all(:,2), [3, fe.sizes.cell]);
    patch(cell_coo_x, cell_coo_y, max(u) * (param / max(param)));
    if verbosity
       fprintf('done.\n'); 
    end
    
    %% Add solution.
    
    if verbosity
       fprintf('Print solution ... '); 
    end
    hold on
        % Display the solution for all DOF (at an adapted triangulation).
        switch fe.order
            case 1
                trisurf(mesh.cell2vtx, x, y, u, 'edgecolor', 'none');
            case 2
                mesh_ref = Mesh.refineMeshUniform(mesh, 1);
                [~, DOF2vtx_map] = ismember(mesh_ref.vertices, ...
                    fe.DOF_maps.DOF_coo, 'rows');
                trisurf(mesh_ref.cell2vtx, ...
                    x(DOF2vtx_map), y(DOF2vtx_map), u(DOF2vtx_map), ...
                    'edgecolor', 'none');
        end
        xlabel('x');
        ylabel('y');
    hold off
    view(3);
    colorbar();
    if verbosity
       fprintf('done.\n'); 
    end
end