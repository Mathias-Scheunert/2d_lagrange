function n = getEdgeNormal(mesh, edge_idx)
    % Calculates the normal vector of a (list of) edge(es).
    %
    % The normal vector is oriented outwards w.r.t. the considered simplex.
    %
    % SYNTAX
    %    n = getEdgeNormal(mesh, edge_idx)
    %
    % INPUT PARAMETER
    %   mesh     ... Struct, containing the mesh information.
    %                For a detailed description of the content of the mesh
    %                struct please read header of Mesh.initMesh.
    %   edge_idx ... Vector [n x 1] of the desired edges indices.
    %
    % OUTPUT PARAMETER
    %   n ... Cell, {n x 2} containing matrix [k x 2]:
    %         global coordinate representation of normal vector(s)
    %         belonging 
    %         first & second columns: x & y component, respectively.

    %% Check imput.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'cell2edg'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isvector(edge_idx), ...
        'edge_idx - vector containing edge indices, expected.');
    
    %% Get cell indices and maps.
    
    % Note that an edge at the interior is related to two simplices.
    [~, cell2_edge_map] = ismember(mesh.cell2edg, edge_idx);
    [edge2cell_idx_map, edge2cell_edge_num] = arrayfun(@(x) ...
                                find(cell2_edge_map == x), ...
                                1:length(edge_idx), ...
                                'UniformOutput', false);
    
    %% Create normal vectors.
          
    % Allocate quantities.
    n_edge = length(edge_idx);
    n = cell(n_edge, 1);
    
    % Loop over edges.
    for ii = 1:n_edge
        
        % Get current simplex number(s).
        cell_num = edge2cell_idx_map{ii};
        
        % Get current simplices coordinates.
        cell_cords = {mesh.cell2cord{cell_num}}.';
        
        % Get edge number within the corresponding simplices.
        edge_num = edge2cell_edge_num{ii};
                
        % Get all edges of current simplices.
        all_cell_edges = arrayfun(@(x) {mesh.cell2edg(x, :)}, ...
                            cell_num);
        
        % Get all edge coordinates for the current simplices.
        all_edges_cords = cellfun(@(x) {reshape({mesh.edge2cord{x(:)}}, ...
                            size(x))}, all_cell_edges);
                        
        % Get current edge coordinates.
        % (Should be equal!)
        edge_vtx_cords = cellfun(@(x, y) {x{y}} , ...
                            all_edges_cords, num2cell(edge_num));
                        
        % Get opposing vertices.
        opposit_vtx = cellfun(@(x, y) {~ismember(x, y, 'rows')}, ...
                         cell_cords, edge_vtx_cords);
        opposit_vtx_cords = cellfun(@(x, y) {x(y, :)}, ...
                               cell_cords, opposit_vtx);
        
        % Calculate normal(s) of current edge and normalize.
        edge_normal = cellfun(@(x) {[-1, 1] .* fliplr(diff(x))}, ...
                         edge_vtx_cords);
        edge_normal = cellfun(@(x) {x / norm(x, 2)}, edge_normal);
        
        % Check their orientation by using another edge of current simplex.
        ref_edge = cellfun(@(x, y) {x(2,:) - y} , ...
                      edge_vtx_cords, opposit_vtx_cords);
        dot_prod = cell2mat(cellfun(@(x, y) {dot(x, y) < 0}, ...
                      edge_normal, ref_edge));
        
        % If required, made them pointing in opposit direction.
        if any(dot_prod)
            edge_normal(dot_prod) = cellfun(@(x) {x .* -1}, ...
                                       edge_normal(dot_prod));
        end

        % Gererate output.
        n{ii} = edge_normal;
    end
end