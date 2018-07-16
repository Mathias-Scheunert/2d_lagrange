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
    %   ymin & ymax: includes the outermost left and right
    %       nodes/vertices/corner points.
    %   xmin & xmax: includes all DOF at the boundary without the corners
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'vertices', 'edge2vtx', 'bnd_edge_xmin'})), ...
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
    idx_xmin_vtx = mesh.edge2vtx(mesh.bnd_edge_xmin, :);
    idx_xmin_vtx = unique(idx_xmin_vtx(:));
    idx_xmax_vtx = mesh.edge2vtx(mesh.bnd_edge_xmax, :);
    idx_xmax_vtx = unique(idx_xmax_vtx(:));
    idx_ymin_vtx = mesh.edge2vtx(mesh.bnd_edge_ymin, :);
    idx_ymin_vtx = unique(idx_ymin_vtx(:));
    idx_ymax_vtx = mesh.edge2vtx(mesh.bnd_edge_ymax, :);
    idx_ymax_vtx = unique(idx_ymax_vtx(:));
    
    % Identify corner point w.r.t vtx identifier.
    vtx_xmin_ymin =  intersect(idx_xmin_vtx, idx_ymin_vtx);
    vtx_xmax_ymin =  intersect(idx_xmax_vtx, idx_ymin_vtx);
    vtx_xmin_ymax =  intersect(idx_xmin_vtx, idx_ymax_vtx);
    vtx_xmax_ymax =  intersect(idx_xmax_vtx, idx_ymax_vtx);
    
    % Remove corner vertices from xmin and xmax bnd.
    remove_idx_min = (idx_xmin_vtx == vtx_xmin_ymin) | ...
        (idx_xmin_vtx == vtx_xmin_ymax);
    remove_idx_max = (idx_xmax_vtx == vtx_xmax_ymin) | ...
        (idx_xmax_vtx == vtx_xmax_ymax);
    idx_xmin_vtx(remove_idx_min) = [];
    idx_xmax_vtx(remove_idx_max) = [];
    
    switch fe.order
        case 1
            % Summarize.
            bnd_DOF_xmin = idx_xmin_vtx;
            bnd_DOF_xmax = idx_xmax_vtx;
            bnd_DOF_ymin = idx_ymin_vtx;
            bnd_DOF_ymax = idx_ymax_vtx;
            [idx_xmin_edg, idx_xmax_edg, idx_ymin_edg, idx_ymax_edg] = ...
                deal([]);

        case 2
            % w.r.t. edges.
            % Get edge index from the already identified bnd edge list and
            % increase the index w.r.t. the number of vertices.
            idx_xmin_edg = n_vtx + find(mesh.bnd_edge_xmin);
            idx_xmax_edg = n_vtx + find(mesh.bnd_edge_xmax);
            idx_ymin_edg = n_vtx + find(mesh.bnd_edge_ymin);
            idx_ymax_edg = n_vtx + find(mesh.bnd_edge_ymax);
            
            % Summarize.
            bnd_DOF_ymin = [idx_ymin_vtx; idx_ymin_edg];
            bnd_DOF_ymax = [idx_ymax_vtx; idx_ymax_edg];
            bnd_DOF_xmin = [idx_xmin_vtx; idx_xmin_edg];
            bnd_DOF_xmax = [idx_xmax_vtx;  idx_xmax_edg];

        otherwise
            error('Unsupported order of Lagrangian elements.');  
    end
    
    % Summarize.
    bnd = struct();
    bnd.bnd_DOF = {bnd_DOF_xmin; bnd_DOF_xmax; bnd_DOF_ymin; bnd_DOF_ymax};
    bnd_DOF = vertcat(bnd.bnd_DOF{:});
    inner_DOF = (1:fe.sizes.DOF).';
    inner_DOF(bnd_DOF) = [];
    bnd.inner_DOF = inner_DOF;
    bnd.n_inner_DOF = length(bnd.inner_DOF);
    bnd.n_bnd_DOF = fe.sizes.DOF - bnd.n_inner_DOF;
        
    %% Map bnd-DOF to their coordinates.
    
    % Collect all DOF positions from vtx or edge midpoint coordinates.    
    bnd.bnd_DOF_coo = {...
    [mesh.vertices(idx_xmin_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_xmin_edg - fe.sizes.vtx), ...
        'UniformOutput', false))];
    [mesh.vertices(idx_xmax_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_xmax_edg - fe.sizes.vtx), ...
        'UniformOutput', false))];
    [mesh.vertices(idx_ymin_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_ymin_edg - fe.sizes.vtx), ...
        'UniformOutput', false))];
    [mesh.vertices(idx_ymax_vtx, :); cell2mat(cellfun(@mean, ...
        mesh.edge2cord(idx_ymax_edg - fe.sizes.vtx), ...
        'UniformOutput', false))]...
        };
end