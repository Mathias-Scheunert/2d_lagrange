function bnd = getBndDOF(fe, mesh)
    % Provides all DOF belonging to boundaries of a rectangular domain.
    %
    % SYNTAX
    %   fe = treatDirichlet(fe, mesh, sol, bnd)
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %
    % OUTPUT PARAMETER
    %   bnd ... Struct, boundary DOF information.
    %
    % REMARKS
    %   Note, that the corner points occur mutiply as each of them is 
    %   shared by two edges!
    %
    %   So the following convention is used:
    %   top & bottom: includes the outermost left and right
    %       nodes/vertices/corner points.
    %   left & rigth: includes all DOF at the boundary without the corners
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'vertices', 'edge2vtx', 'bnd_edge_bot'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(any(strcmp(mesh.type, {'cube', 'rhomb'})), ...
        'mesh.type is not supported. Bnd information has to obtained from external.');
    
    %% Obtain all DOF, belonging to the boundaries.
    
    % Get sizes.
    n_vtx = fe.sizes.vtx;
    
    % w.r.t. vertices.
    % Excluding vertex index from the already identified bnd edge list.
    % Note, that vertex indices occure multiply (as two adjacent domain
    % edges share the corner vertex)!
    idx_bot_vtx = mesh.edge2vtx(mesh.bnd_edge_bot, :);
    idx_bot_vtx = unique(idx_bot_vtx(:));
    idx_top_vtx = mesh.edge2vtx(mesh.bnd_edge_top, :);
    idx_top_vtx = unique(idx_top_vtx(:));
    idx_left_vtx = mesh.edge2vtx(mesh.bnd_edge_left, :);
    idx_left_vtx = unique(idx_left_vtx(:));
    idx_right_vtx = mesh.edge2vtx(mesh.bnd_edge_right, :);
    idx_right_vtx = unique(idx_right_vtx(:));
    
    % Identify corner point w.r.t vtx identifier.
    vtx_lb =  intersect(idx_left_vtx, idx_bot_vtx);
    vtx_rb =  intersect(idx_right_vtx, idx_bot_vtx);
    vtx_lt =  intersect(idx_left_vtx, idx_top_vtx);
    vtx_rt =  intersect(idx_right_vtx, idx_top_vtx);
    
    % Remove corner vertices from left and right bnd.
    remove_idx_left = (idx_left_vtx == vtx_lb) | (idx_left_vtx == vtx_lt);
    remove_idx_right = (idx_right_vtx == vtx_rb) | (idx_right_vtx == vtx_rt);
    idx_left_vtx(remove_idx_left) = [];
    idx_right_vtx(remove_idx_right) = [];
    
    switch fe.order
        case 1
            % Summarize.
            bnd_DOF_bot = idx_bot_vtx;
            bnd_DOF_top = idx_top_vtx;
            bnd_DOF_left = idx_left_vtx;
            bnd_DOF_right = idx_right_vtx; 
            [idx_bot_edg, idx_top_edg, idx_left_edg, idx_right_edg] = ...
                deal([]);

        case 2
            % w.r.t. edges.
            % Get edge index from the already identified bnd edge list and
            % increase the index w.r.t. the number of vertices.
            idx_bot_edg = n_vtx + find(mesh.bnd_edge_bot);
            idx_top_edg = n_vtx + find(mesh.bnd_edge_top);
            idx_left_edg = n_vtx + find(mesh.bnd_edge_left);
            idx_right_edg = n_vtx + find(mesh.bnd_edge_right);
            
            % Summarize.
            bnd_DOF_bot = [idx_bot_vtx; idx_bot_edg];
            bnd_DOF_top = [idx_top_vtx; idx_top_edg];
            bnd_DOF_left = [idx_left_vtx; idx_left_edg];
            bnd_DOF_right = [idx_right_vtx;  idx_right_edg];

        otherwise
            error('Unsupported order of Lagrangian elements.');  
    end
    
    % Summarize.
    bnd = struct();
    bnd.bnd_DOF = {bnd_DOF_bot; bnd_DOF_top; bnd_DOF_left; bnd_DOF_right};
    bnd_DOF = vertcat(bnd.bnd_DOF{:});
    inner_DOF = (1:fe.sizes.DOF).';
    inner_DOF(bnd_DOF) = [];
    bnd.inner_DOF = inner_DOF;
    bnd.n_inner_DOF = length(bnd.inner_DOF);
    bnd.n_bnd_DOF = fe.sizes.DOF - bnd.n_inner_DOF;
        
    %% Map bnd-DOF to their coordinates.
    
    % Collect all DOF positions from vtx or edge midpoint coordinates.    
    bnd.bnd_DOF_coo = {...
    [mesh.vertices(idx_bot_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_bot_edg - fe.sizes.vtx), ...
        'UniformOutput', false))];
    [mesh.vertices(idx_top_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_top_edg - fe.sizes.vtx), ...
        'UniformOutput', false))];
    [mesh.vertices(idx_left_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_left_edg - fe.sizes.vtx), ...
        'UniformOutput', false))];
    [mesh.vertices(idx_right_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_right_edg - fe.sizes.vtx), ...
        'UniformOutput', false))]...
        };
end