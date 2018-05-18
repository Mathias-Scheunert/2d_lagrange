function [] = plotMesh(mesh, params, debug)
    % Visualize the 2d FE mesh.
    % 
    % SYNTAX
    %   [] = plotMesh(mesh[, params, debug])
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %
    % OPTIONAL PARAMETER
    %   params ... Vector, containing the piecewise constant parameter
    %              information, associated with the cells' face.
    %   debug  ... Logical, denoting if additional information will be
    %              visualized.
    
    %% Check input
    
    assert(isstruct(mesh) && isfield(mesh, 'cell2vtx'), ...
        'mesh - struct expected.');
    
    n_cells = length(mesh.cell2vtx);
    if nargin < 3
        debug = false;
        assert(islogical(debug), ...
            'debug - logical expected.');
    end
    if nargin < 2
        params = NaN * zeros(n_cells, 1);
    end
        
    %% Create figure.
    
    figure();
    axis('equal');
    set(gca, 'Ydir', 'reverse')
    xlim([min(mesh.vertices(:,1)), max(mesh.vertices(:,1))]);
    ylim([min(mesh.vertices(:,2)), max(mesh.vertices(:,2))]);
    if all(~isnan(params)) && length(unique(params)) ~= 1
        col_range = [min(params), max(params)] / max(params);
        caxis(col_range);
    else
        caxis('auto');
    end
    if debug
        drawnow();
    end
    
    %% Add simplices.
    
    cell_coo_all = cell2mat(mesh.cell2cord);
    cell_coo_x = reshape(cell_coo_all(:,1), [3, n_cells]);
    cell_coo_y = reshape(cell_coo_all(:,2), [3, n_cells]);
    pause_time = 5/n_cells;
    hold on
    for kk = 1:n_cells
        patch(cell_coo_x(:,kk), cell_coo_y(:,kk), params(kk));
        title(sprintf('num cells: %2.d', kk));
        if debug
            pause(pause_time);
        end
    end
    hold off
    drawnow();
    
    %% Add boundary adges.
    
    if isfield(mesh, 'bnd_edge') && debug
        bnd_edge_cords = mesh.edge2cord(mesh.bnd_edge);
        hold on
        for ii = 1:length(bnd_edge_cords)
            line(bnd_edge_cords{ii}(:,1), bnd_edge_cords{ii}(:,2), ...
                'LineWidth', 4, ...
                'Color', [.8 .8 .8]);
        end
        hold off
    end
    
    %% Add edge orientation.
    
    if isfield(mesh, 'edge2cord') && debug
        hold on
        pos = cell2mat(cellfun(@(x) {x(1,:) + diff(x, 1) / 2}, mesh.edge2cord));
        len = cell2mat(cellfun(@(x) {diff(x, 1)}, mesh.edge2cord));
        scale = 1/3;
        quiver(pos(:,1), pos(:,2), len(:,1), len(:,2), scale, ...
            'Color', 'green', ...
            'LineWidth', 2, ...
            'MaxHeadSize', 1 / norm(len));
        hold off
    end
end
