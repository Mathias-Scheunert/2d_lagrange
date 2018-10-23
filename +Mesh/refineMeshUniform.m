function mesh = refineMeshUniform(mesh, ref_num)
    % Refines mesh by splitting up one triangle into four triangles.
    %
    % This function only acts in the initial mesh information, i.e. the
    % vertices and the cell2vertex list.
    % 
    % SYNTAX
    %   mesh = refineMeshUniform(mesh, ref_num)
    %
    % INPUT PARAMETER
    %   mesh    ... Struct, containing the mesh information.
    %               For a detailed description of the content of the mesh
    %               struct please read header of Mesh.initMesh.
    %   ref_num ... Scalar, denoting the number of refinement steps.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing the refined mesh information.
    %
    % REMARKS
    %   By splitting the triangles, the new nodes, triangles and edges are
    %   placed inbetween the old ones!

    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    
    assert(isscalar(ref_num) && round(ref_num) == ref_num, ...
        'ref_num - integer, denoting number of refinements, expected.');
    
    if ~any(strcmp(mesh.type, {'cube', 'rhomb', 'basic'}))
        error(['Uniform refinement not supported for external meshes. ',...
            'Rather use the functionality of the external mesh generator.']);
        % TODO: Handle splitting up the boundary edge info.
    end
           
    %% Split up triangles.
    
    for i = 1:ref_num
        % Edges are halved and the new points will form the vertices of
        % four new triangles.
        % To make shure that the outermost edge orientation is kept the
        % following principle is used.
        % for the outermost edges:
        % vert.-2-         vert.-2-
        %   edge      new    edge
        %   1', 2'   ->   1, 2, 3
        %   1', 3'   ->   1, 4, 6
        %   2', 3'   ->   3, 5, 6
        % i.e. 
        %   1' -> 1
        %   2' -> 3
        %   3' -> 6
        %
        % vert.-2-        vert.-2-
        %   cell      new   cell
        % 1', 2', 3'  ->  1, 2, 4
        %             ->  2, 3, 5
        %             ->  2, 4, 5
        %             ->  4, 5, 6

        % Get edge midpoints for each cell.   
        edge_mid = cell2mat([...
            cellfun(@(x) ...
            {[x(1,1)/2 + x(2,1)/2; x(1,2)/2 + x(2,2)/2]}, ...
            mesh.cell2cord), ...
            cellfun(@(x) ...
            {[x(1,1)/2 + x(3,1)/2; x(1,2)/2 + x(3,2)/2]}, ...
            mesh.cell2cord), ...
            cellfun(@(x) ...
            {[x(2,1)/2 + x(3,1)/2; x(2,2)/2 + x(3,2)/2]}, ...
            mesh.cell2cord) ...     
            ]);

        % Reshape expression to coincide with mesh.cell2cord structure.
        n_cell_old = length(edge_mid) / 2;
        edge_mid = mat2cell(edge_mid, 2 + zeros(n_cell_old, 1), 3);
        edge_mid = cellfun(@transpose, edge_mid, 'UniformOutput', false);

        % Redefine cell2cord list by incorporating the new cells.
        % That is, one former cell will be replaced by four new cells.
        % Therefore, use the above mentioned vtx2cell ordering.
        cell2cord = cellfun(@(x,y) {...
            [[x(1,1), x(1,2)]; [y(1,1), y(1,2)];[y(2,1), y(2,2)]]; ...
            [[y(1,1), y(1,2)]; [x(2,1), x(2,2)];[y(3,1), y(3,2)]]; ...
            [[y(1,1), y(1,2)]; [y(2,1), y(2,2)];[y(3,1), y(3,2)]];
            [[y(2,1), y(2,2)]; [y(3,1), y(3,2)];[x(3,1), x(3,2)]]}, ...
            mesh.cell2cord, edge_mid, 'UniformOutput', false);

        % Define 'local' vertex lists for each cell.
        vert_list_local = cellfun(@(x) ...
            [x{1}(1,1), x{1}(1,2); 
             x{1}(2,1), x{1}(2,2);
             x{2}(2,1), x{2}(2,2);
             x{3}(2,1), x{3}(2,2);
             x{3}(3,1), x{3}(3,2);
             x{4}(3,1), x{4}(3,2)], cell2cord, 'UniformOutput', false);

        % Generate 'global' vertex list.
        vert_list_local = cell2mat(vert_list_local);
        [vert_list, ~, ind_glob2loc] = unique(vert_list_local, ...
            'rows', 'stable');

        % Define 'local' cell list.
        cell_list_local = [1, 2, 4;
                           2, 3, 5;
                           2, 4, 5;
                           4, 5, 6];
        cell_list_local_tmp = cell_list_local.';
        cell_list_local_tmp = cell_list_local_tmp(:);

        % Reshape indGlob2Loc to fit the mesh.cell2cord structure.
        ind_glob2loc_tmp = mat2cell(ind_glob2loc, 6 + zeros(n_cell_old, 1), 1);

        % Generate 'global' cell list.
        % (formal lexicographic order is discarted)
        cell_list = cellfun(@(x) ...
            reshape(x(cell_list_local_tmp), 3, 4).', ...
            ind_glob2loc_tmp, 'UniformOutput', false);
        cell_list = cell2mat(cell_list);
        
        % Make sure, that cell2vtx is sorted ascendingly.
        cell_list = sort(cell_list, 2);

        if strcmp(mesh.type, 'basic')
            % Initialize new (avoid mixing of refined and unrefined info).
            mesh = struct();
            mesh.type = 'basic';
            mesh.vertices = vert_list;
            mesh.cell2vtx = cell_list;
            
        else
            % Override existing information.
            mesh.vertices = vert_list;
            mesh.cell2vtx = cell_list;
        end
        
        if ~strcmp(mesh.type, 'basic')
            % Expand parameter domain vector by inserting new cells.
            % (Relying on the 1 -> 4 rule [see cell2cord definition])
            mesh.parameter_domain = kron(mesh.parameter_domain(:), ...
                                        ones(4, 1));
        
            % Update coordinate relations.
            mesh = Mesh.appendEdgeInfo(mesh);
            mesh = Mesh.appendCoordInfo(mesh);
            mesh = Mesh.appendBndInfo(mesh);
        end
    end
end
