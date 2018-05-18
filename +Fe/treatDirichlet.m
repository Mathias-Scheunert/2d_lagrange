function sol = treatDirichlet(fe, mesh, sol, bnd)
    % Adapts the FE linear system to handle Dirichlet boundary conditions.
    %
    % SYNTAX
    %   fe = treatDirichlet(fe, mesh, sol, bnd)
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %   sol  ... Struct, containing the information of the current physical
    %            problem to be solved numerically, i.e. rhs vector, system 
    %            matrix, interpolation operator.
    %   bnd  ... Struct, containing the boundary condition information.
    %
    % OUTPUT PARAMETER
    %   sol ... Struct, adapted sol struct (see INPUT PARAMETER).
    %
    % REMARKS
    %   The order of dirichlet bnd values given in bnd.val (ONLY in case of
    %   length(bnd.val) == 4) is:
    %   bnd.val = [bot; top; left; right]
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(sol) && all(isfield(sol, { 'A', 'b'})), ...
        'sol - struct, containing info about linear system to be solved, expected.');
    assert(isstruct(bnd) && all(isfield(bnd, {'type', 'val'})), ...
        'bnd - struct, containing the boundary condition information, expected.');
    assert(strcmp(bnd.type, 'dirichlet'), ...
        'bnd.type - only handling of Dirichlet values are supported here.');
 
    %% Obtain all DOF, belonging to the boundaries.
    
    % Get sizes.
    n_vtx = fe.sizes.vtx;
    
    % w.r.t. vertices.
    idx_bot_vtx = mesh.edge2vtx(mesh.bnd_bot, :);
    idx_bot_vtx = unique(idx_bot_vtx(:));
    idx_top_vtx = mesh.edge2vtx(mesh.bnd_top, :);
    idx_top_vtx = unique(idx_top_vtx(:));
    idx_left_vtx = mesh.edge2vtx(mesh.bnd_left, :);
    idx_left_vtx = unique(idx_left_vtx(:));
    idx_right_vtx = mesh.edge2vtx(mesh.bnd_right, :);
    idx_right_vtx = unique(idx_right_vtx(:));
    
    switch fe.order
        case 1
            % Summarize.
            bnd_DOF_bot = idx_bot_vtx;
            bnd_DOF_top = idx_top_vtx;
            bnd_DOF_left = idx_left_vtx;
            bnd_DOF_right = idx_right_vtx;  
            
        case 2
            % w.r.t. edges.
            idx_bot_edg = find(mesh.bnd_bot);
            idx_top_edg = find(mesh.bnd_top);
            idx_left_edg = find(mesh.bnd_left);
            idx_right_edg = find(mesh.bnd_right);
            
            % Summarize.
            bnd_DOF_bot = [idx_bot_vtx; n_vtx + idx_bot_edg];
            bnd_DOF_top = [idx_top_vtx; n_vtx + idx_top_edg];
            bnd_DOF_left = [idx_left_vtx; n_vtx + idx_left_edg];
            bnd_DOF_right = [idx_right_vtx; n_vtx + idx_right_edg];

        otherwise
            error('Unsupported order of Lagrangian elements.');  
    end
    
    % Summarize.
    bnd_DOF = [bnd_DOF_bot; bnd_DOF_top; bnd_DOF_left; bnd_DOF_right];
    inner_DOF = (1:fe.sizes.DOF).';
    inner_DOF(bnd_DOF) = [];
    
    % Check consistency for given bnd constraint.
    if isvector(bnd.val) && length(bnd.val) == 4
        assert(all(bnd.val(1) == bnd.val(2:end)), ...
            ['For uniform Dirichlet bnd-cond all values must be equal. ', ...
            'Otherwise inconsistencies occur at the corners of the grid.']);
    else
        error('Unsupported Dirichlet boundary value distribution.');
        % TODO: implement.
    end
    
    %% Reduce the linear system.
       
    % RHS vector.
    if any(bnd.val ~= 0)
        % To couple the inhomogeneouse Dirichlet values with the other DOF
        % the rhs needs to be modyfied.
        b_dirichlet = zeros(fe.sizes.DOF ,1);
        b_dirichlet(bnd_DOF_bot) = bnd.val(1);
        b_dirichlet(bnd_DOF_top) = bnd.val(2);
        b_dirichlet(bnd_DOF_left) = bnd.val(3);
        b_dirichlet(bnd_DOF_right) = bnd.val(4);
        b_dirichlet = sol.A * b_dirichlet;
        sol.b = sol.b - b_dirichlet;
    end
    sol.b = sol.b(inner_DOF);
    
    % System matrix.
    sol.A = sol.A(inner_DOF, inner_DOF);
    
    %% Append Dirichlet bnd info.
    
    sol.dirichlet.bnd_DOF = bnd_DOF;
    sol.dirichlet.inner_DOF = inner_DOF;
    sol.dirichlet.bnd_DOF_bot = bnd_DOF_bot;
    sol.dirichlet.bnd_DOF_top = bnd_DOF_top;
    sol.dirichlet.bnd_DOF_left = bnd_DOF_left;
    sol.dirichlet.bnd_DOF_right = bnd_DOF_right;
    sol.dirichlet.val = bnd.val;
end