function bnd = getBndDOF(fe, mesh, bnd)
    % Provides all DOF belonging to boundaries of the given 2D domain.
    %
    % SYNTAX
    %   bnd = getBndDOF(fe, mesh, bnd)
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %   bnd  ... Struct, containing the preset of the BC types.
    %
    % OUTPUT PARAMETER
    %   bnd ... bnd Struct appanded by the boundary DOF information.
    %
    % REMARKS
    %   Note, that a corner point of adjacent boundary parts occur mutiply
    %   as it is shared by two edges.

    %% Check input.

    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - Struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'vertices', 'edge2vtx', 'bnd_edge'})), ...
        'mesh - Appended struct, including edge and mapping information, expected.');
    assert(any(strcmp(mesh.type, {'cube', 'disc', 'gmsh_create', 'gmsh_load'})), ...
        'mesh.type is not supported.');
    assert(isstruct(bnd) && all(isfield(bnd, {'name'})), ...
        'bnd - Struct, including general information of bnd types expected.');

    %% Check input.

    % Find similar entries within both lists.
    match_bnd = ismember(bnd.name, mesh.bnd_edge_part_name);
    match_mesh = ismember(mesh.bnd_edge_part_name, bnd.name);

    % Check if entries are missing.
    if ~all(match_bnd)
        warning(['There are boundary part names stored in the bnd ', ...
        'struct which are not given by the mesh struct. Please make ', ...
        'sure that all boundary parts of the mesh are related to an ', ...
        'appropriate boundary condition while setting up the ', ...
        'DRIVE_ script.\n', ...
        'Following boundaries of "bnd" will be removed:', ...
        repmat('\n\t"%s"', 1, length(find(~match_bnd)))], ...
        bnd.name{~match_bnd});

        % Remove obsolete entries.
        bnd.name = bnd.name(match_bnd);
        bnd.val = cellfun(@(x) {x(match_bnd)}, bnd.val);
    end

    if ~all(match_mesh)
        warning(['There are boundary names stored in the mesh ', ...
        'which are not given in the bnd struct. Please make sure that ', ...
        'all boundary parts of the mesh are related to an appropriate ', ...
        'boundary condition while setting up the DRIVE_ script.\n', ...
        'Following boundaries of "mesh" will be ignored:', ...
        repmat('\n\t"%s"', 1, length(find(~match_mesh)))], ...
        mesh.bnd_edge_part_name{~match_mesh});

        % Remove obsolete entries.
        bnd_edge2idx = mesh.bnd_edge_part == 1:length(mesh.bnd_edge_part_name);
        mesh.bnd_edge_part_name = mesh.bnd_edge_part_name(match_mesh);
        mesh.bnd_edge = bnd_edge2idx * double(match_mesh);
        mesh.bnd_edge_part = bnd_edge2idx(:,match_mesh) * ...
                                (1:length(find(match_mesh))).';
    end

    % Sort remaining information of bnd to match the order of mesh.
    [~, sort_like_in_mesh] = ismember(mesh.bnd_edge_part_name, bnd.name);
    bnd.name = bnd.name(sort_like_in_mesh);
    bnd.val = cellfun(@(x) {x(sort_like_in_mesh)}, bnd.val);

    %% Obtain all DOF, belonging to the boundaries.

    % Get sizes.
    n_vtx = fe.sizes.vtx;
    n_bnd_part = length(mesh.bnd_edge_part_name);

    % w.r.t. vertices.
    % Excluding vertex index from the already identified bnd edge list.
    % Note, that vertex indices occure multiply (as two adjacent domain
    % edges share the corner vertex)!
    bnd_part_vtx_idx = cell(n_bnd_part, 1);
    for ii = 1:n_bnd_part
        bnd_part_vtx_idx{ii} = mesh.edge2vtx(mesh.bnd_edge_part == ii,:);
        bnd_part_vtx_idx{ii} = unique(bnd_part_vtx_idx{ii}(:));
    end

    switch fe.order
        case 1
            % Summarize.
            bnd_DOF_list = bnd_part_vtx_idx;
            bnd_part_edge_idx = cell(n_bnd_part, 1);

        case 2
            % w.r.t. edges.
            % Get edge index from the already identified bnd edge list and
            % increase the index w.r.t. the number of vertices.
            bnd_part_edge_idx = cell(n_bnd_part, 1);
            for ii = 1:n_bnd_part
                bnd_part_edge_idx{ii} = n_vtx + ...
                                            find(mesh.bnd_edge_part == ii);
            end

            % Summarize.
            bnd_DOF_list = cellfun(@(x, y) {[x; y]}, bnd_part_vtx_idx, ...
                              bnd_part_edge_idx);

        otherwise
            error('Unsupported order of Lagrangian elements.');
    end

    % Summarize.
    bnd_DOF = struct();
    bnd_DOF.bnd_DOF = bnd_DOF_list;
    tmp_DOF = vertcat(bnd_DOF.bnd_DOF{:});
    inner_DOF = (1:fe.sizes.DOF).';
    inner_DOF(tmp_DOF) = [];
    bnd_DOF.inner_DOF = inner_DOF;
    bnd_DOF.n_inner_DOF = length(bnd_DOF.inner_DOF);
    bnd_DOF.n_bnd_DOF = fe.sizes.DOF - bnd_DOF.n_inner_DOF;

    %% Map bnd-DOF to their coordinates.

    % Collect all DOF positions from vtx or edge midpoint coordinates.
    bnd_DOF_coo = cell(n_bnd_part, 1);
    for ii = 1:n_bnd_part
        bnd_DOF_coo{ii} = [mesh.vertices(bnd_part_vtx_idx{ii}, :); ...
                            cell2mat(cellfun(@mean, ...
                            mesh.edge2cord(bnd_part_edge_idx{ii} - fe.sizes.vtx), ...
                            'UniformOutput', false))];
    end
    bnd_DOF.bnd_DOF_coo = bnd_DOF_coo;

    %% Appand bnd struct.

    bnd.bndDOF = bnd_DOF;
end
