function bnd = getBndDOF(fe, mesh)
    % Provides all DOF belonging to boundaries.
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
    %   See also: Fe.treatDirichlet, Fe.treatNeumann
    %
    % TODO: it may be wise to fix that multiply occurence!
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    
    %% Obtain all DOF, belonging to the boundaries.
    
    % Get sizes.
    n_vtx = fe.sizes.vtx;
    
    % w.r.t. vertices.
    % Excluding vertex index from the already identified bnd edge list.
    idx_bot_vtx = mesh.edge2vtx(mesh.bnd_edge_bot, :);
    idx_bot_vtx = unique(idx_bot_vtx(:));
    idx_top_vtx = mesh.edge2vtx(mesh.bnd_edge_top, :);
    idx_top_vtx = unique(idx_top_vtx(:));
    idx_left_vtx = mesh.edge2vtx(mesh.bnd_edge_left, :);
    idx_left_vtx = unique(idx_left_vtx(:));
    idx_right_vtx = mesh.edge2vtx(mesh.bnd_edge_right, :);
    idx_right_vtx = unique(idx_right_vtx(:));
    % Note, that vertex indices occure multiply (as two adjacent domain
    % edges share the corner vertex)!
    n_bnd_vtx_full = length(idx_bot_vtx) + length(idx_top_vtx) + ...
        length(idx_left_vtx) + length(idx_right_vtx);
    
    switch fe.order
        case 1
            % Summarize.
            bnd_DOF_bot = idx_bot_vtx;
            bnd_DOF_top = idx_top_vtx;
            bnd_DOF_left = idx_left_vtx;
            bnd_DOF_right = idx_right_vtx;  
            
        case 2
            % w.r.t. edges.
            % Get edge index from the already identified bnd edge list.
            idx_bot_edg = find(mesh.bnd_edge_bot);
            idx_top_edg = find(mesh.bnd_edge_top);
            idx_left_edg = find(mesh.bnd_edge_left);
            idx_right_edg = find(mesh.bnd_edge_right);
            
            % Summarize.
            bnd_DOF_bot = [idx_bot_vtx; n_vtx + idx_bot_edg];
            bnd_DOF_top = [idx_top_vtx; n_vtx + idx_top_edg];
            bnd_DOF_left = [idx_left_vtx; n_vtx + idx_left_edg];
            bnd_DOF_right = [idx_right_vtx; n_vtx + idx_right_edg];

        otherwise
            error('Unsupported order of Lagrangian elements.');  
    end
    
    % Summarize.
    bnd = struct();
    bnd.DOF_bot = bnd_DOF_bot; 
    bnd.DOF_top = bnd_DOF_top; 
    bnd.DOF_left = bnd_DOF_left;
    bnd.DOF_right = bnd_DOF_right;
    bnd.bnd_DOF = [bnd_DOF_bot; bnd_DOF_top; bnd_DOF_left; bnd_DOF_right];
    bnd.n_bnd_DOF_unique = length(unique(bnd.bnd_DOF));
    inner_DOF = (1:fe.sizes.DOF).';
    inner_DOF(bnd.bnd_DOF) = [];
    bnd.inner_DOF = inner_DOF;
    
    %% Add corner point w.r.t vtx identifier.
    
    bnd.vtx_lb =  intersect(idx_left_vtx, idx_bot_vtx);
    bnd.vtx_rb =  intersect(idx_right_vtx, idx_bot_vtx);
    bnd.vtx_lt =  intersect(idx_left_vtx, idx_top_vtx);
    bnd.vtx_rt =  intersect(idx_right_vtx, idx_top_vtx);
    
    %% Map bnd-DOF to their coordinates.
    
    % Collect all DOF positions from vtx or edge midpoint coordinates.
    n_bnd_DOF_full = length(bnd.bnd_DOF);
    bnd_vtx = bnd.bnd_DOF <= fe.sizes.vtx;
    bnd_edge = bnd.bnd_DOF > fe.sizes.vtx;
    bnd_DOF_coo = zeros(n_bnd_DOF_full, 2);
    bnd_DOF_coo(bnd_vtx,:) = mesh.vertices(bnd.bnd_DOF(bnd_vtx),:);
    bnd_DOF_coo(bnd_edge,:) = cell2mat(cellfun(@mean, ...
        mesh.edge2cord(bnd.bnd_DOF(bnd_edge) - fe.sizes.vtx), ...
        'UniformOutput', false));
        
    % Summarize.
    bnd.bnd_DOF_coo = bnd_DOF_coo;
end