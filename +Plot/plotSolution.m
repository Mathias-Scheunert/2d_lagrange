function [] = plotSolution(fe, mesh, u, verbosity)
    % Visualize the solution of the 2D FE fwd problem.
    % 
    % SYNTAX
    %   [] = plotSolution(fe, mesh, u, verbosity)
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
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
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
    if verbosity
       fprintf('Create figure ... '); 
    end
    figure();
    axis('equal');
    set(gca, 'Ydir', 'reverse')
    xlim([min(mesh.vertices(:,1)), max(mesh.vertices(:,1))]);
    ylim([min(mesh.vertices(:,2)), max(mesh.vertices(:,2))]);
    if verbosity
       fprintf('Print mesh ... '); 
    end
    trimesh(mesh.cell2vtx, x(1:fe.sizes.vtx), y(1:fe.sizes.vtx), y(1:fe.sizes.vtx) * 0);
    
    %% Add solution.
    
    if verbosity
       fprintf('Print solution ... '); 
    end
    hold on
        % Display the solution for all DOF at an adapted triangulation.
        tri = delaunay(x, y);
        trisurf(tri, x, y, u, ...
            'edgecolor', 'none');
    hold off
    view(3);
    colorbar();
    if verbosity
       fprintf('done.\n'); 
    end
end