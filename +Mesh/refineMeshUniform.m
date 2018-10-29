function mesh = refineMeshUniform(mesh, ref_num)
    % Refines mesh by splitting up one triangle into four triangles.
    %
    % This function only acts in the initial mesh information, i.e. the
    % vertices, the cell2vtx list and the parameter_domain list.
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
    %   Edges are halved and the new points will form the vertices of
    %   four new triangles.
    %   To make sure that the outermost edge orientation is kept the
    %   following principle is used:
    %       .' entity on the underlying (global) mesh
    %       .  entity on the refined (local) mesh
    %
    %   global relations:
    %    edge2vtx'     edge2vtx
    %      1', 2'  ->   1, 2, 3
    %      2', 3'  ->   3, 5, 6
    %      1', 3'  ->   1, 4, 6
    %
    %     vtx'   vtx
    %      1' ->  1
    %      2' ->  3
    %      3' ->  6
    %
    %      cell2vtx'      cell2vtx
    %      1', 2', 3'  ->  1, 2, 4
    %                  ->  2, 3, 5
    %                  ->  2, 4, 5
    %                  ->  4, 5, 6
    %
    %   local relations:
    %      cell'(edge2vtx')      cell(edge2vtx)
    %         1'(1',2')    ->    1(1,2), 2(1,2)
    %         1'(2',3')    ->    2(2,3), 4(2,3)
    %         1'(1',3')    ->    1(1,3), 4(1,3)

    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    
    assert(isscalar(ref_num) && round(ref_num) == ref_num, ...
        'ref_num - integer, denoting number of refinements, expected.');
           
    %% Split up triangles.
    
    for i = 1:ref_num
        % Get edge midpoints for each cell.   
        cell_edge_mid = cell2mat([...
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
        n_cell_old = length(cell_edge_mid) / 2;
        cell_edge_mid = mat2cell(cell_edge_mid, 2 + zeros(n_cell_old, 1), 3);
        cell_edge_mid = cellfun(@transpose, cell_edge_mid, ...
                                'UniformOutput', false);

        % Redefine cell2cord list by incorporating the new cells.
        % That is, one former cell will be replaced by four new cells.
        % Therefore, use the above mentioned vtx2cell ordering.
        cell2cord = cellfun(@(x,y) {...
            [[x(1,1), x(1,2)]; [y(1,1), y(1,2)];[y(2,1), y(2,2)]]; ...
            [[y(1,1), y(1,2)]; [x(2,1), x(2,2)];[y(3,1), y(3,2)]]; ...
            [[y(1,1), y(1,2)]; [y(2,1), y(2,2)];[y(3,1), y(3,2)]];
            [[y(2,1), y(2,2)]; [y(3,1), y(3,2)];[x(3,1), x(3,2)]]}, ...
            mesh.cell2cord, cell_edge_mid, 'UniformOutput', false);

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

        % Reshape ind_glob2loc to fit the mesh.cell2cord structure.
        ind_glob2loc_tmp = mat2cell(ind_glob2loc, 6 + zeros(n_cell_old, 1), 1);

        % Generate 'global' cell list.
        % (formal lexicographic order is discarted)
        cell_list = cellfun(@(x) ...
                        reshape(x(cell_list_local_tmp), 3, 4).', ...
                        ind_glob2loc_tmp, 'UniformOutput', false);
        cell_list = cell2mat(cell_list);
        
        % Generate local to global cell map.
        cell_loc2glob = kron((1:n_cell_old).', ...
                            ones(size(cell_list_local, 1), 1));
                        
        % Generate local to global edge map.
        edge_cell_loc2glo = {[1, 2];
                             [2, 4];
                             [1, 4]};
        % (Note ordering, see Mesh.appendEdgeInfo and file header).
        edge_cell_edge_loc2glo = {[1, 2], [1, 2];
                                  [2, 3], [2, 3];
                                  [1, 3], [1, 3]};
        edge_loc2glo = [edge_cell_loc2glo, edge_cell_edge_loc2glo];
        
        % Make sure, that cell2vtx is sorted ascendingly.
        cell_list = sort(cell_list, 2);
        
        if strcmp(mesh.type, 'basic')
            % Initialize new (avoid mixing of refined and unrefined info).
            mesh = struct();
            mesh.type = 'basic';
            mesh.vertices = vert_list;
            mesh.cell2vtx = cell_list;
            
        else            
            % Add mapping.
            if isfield(mesh, 'cell_loc2glob')
                mesh.cell_loc2glob = [mesh.cell_loc2glob, cell_loc2glob];
                mesh.edge_loc2glo = [mesh.edge_loc2glo; edge_loc2glo];
            else
                mesh.cell_loc2glob = {cell_loc2glob};
                % TODO: check if this info is really required to be stored.
                mesh.edge_loc2glo = {edge_loc2glo};
            end

            % Expand parameter domain vector by inserting new cells.
            mesh.parameter_domain = kron(mesh.parameter_domain(:), ...
                                        ones(size(cell_list_local, 1), 1));
                                    
            % Handle external meshes (i.e. known boundary info).
            if any(strcmp(mesh.type, {'gmsh_create', 'gmsh_load'}))
                
                % Expand essential boundary edge info (which can't be
                % obtained otherwise).
                % Get edge midpoints for each bnd edge.   
                bnd_edge_mid = mesh.edge2cord(mesh.bnd_edge);
                bnd_edge_mid = cell2mat(...
                    cellfun(@(x) {[(x(1,1) + x(2,1))/2, ...
                                   (x(1,2) + x(2,2))/2]}, ...
                                 bnd_edge_mid) ...
                               );
                           
                % Get vtx (w.r.t. refined mesh) id for those points.
                % TODO: avoid unnessecary workload here! May use info from
                %       cell_edge_mid list instead of the following?
                [~, bnd_edge_mid_vtx_locinglo] = ismember(vert_list, ...
                                                     bnd_edge_mid, 'rows');
                [glo_vtx_idx, ~, glo2loc_map] = find(bnd_edge_mid_vtx_locinglo);
                [~, glo2loc_idx] = sort(glo2loc_map);
                bnd_edge_mid_vtx_locinglo = glo_vtx_idx(glo2loc_idx);
                assert(length(bnd_edge_mid_vtx_locinglo) == size(bnd_edge_mid, 1), ...
                    ['Not all boundary edge midpoints could be found in ', ...
                    'vertices list.']);
                               
                % Get edge2vtx list for the unrefined mesh.
                bnd_edge2vtx = mesh.edge2vtx(mesh.bnd_edge,:);
                
                % Get the cooresponding vertices in the refined mesh.
                [bnd_edge_vtx, ~, lin_idx] = unique(bnd_edge2vtx);
                bnd_edge_vtx_cords = mesh.vertices(bnd_edge_vtx, :);
                [~, bnd_edge_vtx_loc2glo] = ismember(bnd_edge_vtx_cords, ...
                                                     vert_list, 'rows');
                
                % Replace vertices index with their number w.r.t. refined
                % mesh.
                bnd_edge2vtx = reshape(bnd_edge_vtx_loc2glo(lin_idx), ...
                                       size(bnd_edge2vtx));
                
                % Include midpoint twice and split up bnd_edge2vtx list.
                bnd_edge2vtx = [bnd_edge2vtx(:, 1) , bnd_edge_mid_vtx_locinglo, ...
                                bnd_edge_mid_vtx_locinglo, bnd_edge2vtx(:, 2)];
                bnd_edge2vtx = reshape(bnd_edge2vtx.', ...
                                   2, size(bnd_edge2vtx, 1) * 2).';
                mesh.bnd_edge2vtx = bnd_edge2vtx;
                
                % Reduce vector to the size of the known bnd edges as it
                % will be extended appropriately within
                % Mesh.appendEdgeInfo.
                % (keep order, as it fits to the bnd_edge order)
                [~, ~, mesh.bnd_edge_part] = find(mesh.bnd_edge_part);
                mesh.bnd_edge_part = kron(mesh.bnd_edge_part, [1; 1]);
            end
            
            % Override existing information.
            mesh.vertices = vert_list;
            mesh.cell2vtx = cell_list;
            
            % Remove unrefined info from mesh (which was't be replaced).
            mesh = rmfield(mesh, {'cell2cord', 'edge2cord', ...
                                  'cell2edg', 'edge2vtx', ...
                                  'bnd_edge'});
        
            % Update coordinates and further mesh relations.
            mesh = Mesh.appendEdgeInfo(mesh);
            mesh = Mesh.appendCoordInfo(mesh);
            mesh = Mesh.appendBndInfo(mesh);
        end
    end
end
