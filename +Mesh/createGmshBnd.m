function [mesh] = createGmshBnd(bnd_in,TX,RX,topo,dom_name,args)
    % Stores mesh information of a mesh created by Gmsh.
    %
    % Supported (tested) Gmsh version: 3.x, MSH file format version 2
    %
    % Therefore,
    %   1) a .geo input file is created
    %   2) a .msh output file is obtained from running Gmsh
    %   3) information from .msh is converted
    %   4) .geo and .msh files are deleted
    %
    % SYNTAX
    %   mesh = createGmsh(bnd, args)
    %
    % INPUT PARAMETER
    %   bnd  ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   args ... Struct, may containing further vectorized (TX, RX, topo)
    %            and scalar (ref) information.
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing vertex2coordinates (vertices),
    %            simplex2vertex (cell2vtx), boundary egde2coordinates
    %            (bnd_edge2vtx), boundary egde 2 model domain boundary
    %            (bnd_edge_...) and simplex parameter domains
    %            (parameter_domains).
    % field parameters: lc_min = 0.05; lc_max = dom_wdt/2; DistMax = dom_wdt; DistMin = 0.8;  Mesh.Algorithm=1


    %% Check Gmsh version.

    gmsh_path = dir('**/gmsh.exe');
    if isempty(gmsh_path)
        error('Gmsh executable could not be found in path.');
    end

    %% Create Gmsh input file.

    gmsh_file = 'homo_mesh_bnd';
    createGmshInput([gmsh_file, '.geo'], bnd_in, TX, RX, topo, dom_name,args.verbosity);

    %% Run Gmsh.
     system([gmsh_path.folder, '/gmsh -2 ', ...
            gmsh_file, '.geo -v 0 -format msh2 ']);

    %% Import mesh information from .msh file.

    mesh = Mesh.loadGmsh([gmsh_file, '.msh'], ...
        'verbosity', args.verbosity, ...
        'ref', args.ref);

    % Override type.
    % TODO: may remove this type completely.
    mesh.type = 'gmsh_create';

    %% Clean up.

    %delete([gmsh_file, '.geo'], [gmsh_file, '.msh']);
end

function createGmshInput(name, bnd_in, TX, RX, topo, dom_name, verbosity)
    % Creates an .geo input file for Gmesh
    %
    % SYNTAX
    %   createGmshInput(name, bnd, TX, RX, topo, dom_name)
    %
    % INPUT PARAMETER
    %   name  ... Char, denoting the file name to be created.
    %   bnd   ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   TX    ... Vector [n x 2], denoting the source position(s).
    %   RX    ... Vector [m x 2], denoting the receiver position(s).
    %   topo  ... Vector [k x 2], denoting the descrete topography.
    %   dom_name ... Char, denoting the physical domain name.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed
    %
    % REMARKS
    %  Note that the 2D grid will be referred to gmesh x and z coordinate:
    %  0,0 ------>
    %      |      x = gmesh-x
    %      |
    %      |
    %      v
    %       y = gmesh-z
    %
    %  Note that gmsh can't handle point which lie on the boundary!
    %  I.g. if a TX/RX point is embedded in the surface but lies somewhere
    %  on a bounding line of that surface, gmsh will not be able to find a
    %  mesh for that situation.
    %  This might be due to the behavior, that gmsh meshes starting from
    %  lines to surfaces to volumes an can't handle additional points in
    %  surface which also occur on the already meshed lines.
    %
    %  To fix that issue, TX/RX points wich are not horizontaly aligned
    %  with a topography point, are expected to lie ON the boundary, such
    %  that the domain boundary will be adapted.
    %% Check input.
    %in initMesh
    if verbosity
       fprintf('Define Gmsh input ... ');
    end

    %% Set domain boundaries.

    TXRX = unique([TX; RX], 'rows');
    if ~isempty(TXRX)

        % Get specific extentions.
        %min_x = min(TXRX(:,1));
        %max_x = max(TXRX(:,1));
        %offset = max_x-min_x;

        % Adapt.
        bnd = [bnd_in(1), ...
               bnd_in(2), ...
               bnd_in(3), ...
               bnd_in(4)];  % +offset all uncommented bc improved bnd
    end

    % Transform boundary coordinates.
    X = [bnd(1:2), fliplr(bnd(1:2))];
    Y = [bnd(3:4); bnd(3:4)];
    Y = Y(:).';

    %% Set up domain basic geometry entities (points, lines).

    % Set maximal element size.
    %                  geometry based                  fix
    % removed and replaced with distance fields
    % h_domain = pick(1, 5 * round(diff(bnd(1:2)) / 25), 10);

    % Create point input ([n x 5] matrix) structure for domain boundaries.
    n_domain = length(X);
    point_domain = zeros(n_domain, 4);
    for i = 1:n_domain
       point_domain(i, 1) = i ;       % point number
       point_domain(i, 2) = X(i);     % gmesh-x cooridnate.
       point_domain(i, 3) = 0;        % gmesh-y cooridnate.
       point_domain(i, 4) = Y(i);     % gmesh-z cooridnate.
       %point_domain(i, 5) = h_domain; % element size around current point
    end

    % Add topography points to point list.
%     if ~isempty(topo)
%
%         % Set maximal element size.
%         % Use min point offset in x direction (~horizontal) & avoid zero.
%         h_topo = max(min(diff(topo(:,1))), 0.1);
%
%         % Create point input ([n x 5] matrix) structure for topography.
%         n_topo = size(topo, 1);
%         n_domain = n_topo + n_domain;
%         point_topo = [(1:n_topo).', topo(:,1), zeros(n_topo, 1), ...
%             topo(:,2), h_topo + zeros(n_topo, 1)];
%
%         % Insert points into domain point list.
%         point_domain = [point_domain(1, :); ...
%                         point_topo; ...
%                         point_domain(2:end, :)];
%         point_domain(:,1) = 1:n_domain;
%     end

    % Add TX/RX positions.
    if ~isempty(TXRX)

        tol_find_TX_in_pd = pick(1, 1e-3, 1e-6);
        [~, idx_txrx_in_pd] = ismembertol(...
                                       point_domain(:, [2, 4]), TXRX, ...
                                       tol_find_TX_in_pd, 'ByRows', true);
        %point_domain(txrx_in_pd, 5) = h_TXRX;
        idx_txrx_not_in_pd = ~ismember(...
            1:size(TXRX, 1), idx_txrx_in_pd(idx_txrx_in_pd ~= 0));

        % Handle all TX/RX points that are not included in topo.
        if any(idx_txrx_not_in_pd)
            n_TXRX = length(find(idx_txrx_not_in_pd));
            % Check if left TX/RX points are at least horizontally aligned
            % with topography information.
            if ~isempty(topo) && any(ismember(...
                    TXRX(idx_txrx_not_in_pd, 1), ...
                    topo(:, 1)))
                % Get corresponding point index in topo.
                [idx_in_topo, ~] = find(ismember(...
                                        topo(:, 1), ...
                                        TXRX(idx_txrx_not_in_pd, 1) ...
                                        ));
                % Get corresponding point index in TXRX.
                [idx_in_TXRX, ~] = find(ismember(...
                                        TXRX(idx_txrx_not_in_pd, 1), ...
                                        topo(:, 1) ...
                                        ));
                % Check, if any TX/RX point lies above the surface.
                if any(TXRX(idx_in_TXRX, 2) < topo(idx_in_topo, 2))
                    error(['TX/RX-Point lies outside the model ', ...
                        'domain which is specified by the topography.']);
                end

                % Otherwise add points to interior of the model domain.
                % (finer mesh is needed)
                n_TXRX_intern = length(idx_in_topo);
                n_TXRX_add = n_TXRX - n_TXRX_intern;
                tmp_idx = find(idx_txrx_not_in_pd);
                idx_internal = tmp_idx(idx_in_TXRX);
                % Make sure that the ongoing numbering for the internal
                % points starts at the end of all points in point_domain
                % list (also including those TX/RX which are not handled
                % yet)
                n_domain = n_domain + n_TXRX_add;
                point_internal = [n_domain + (1:n_TXRX_intern).', ...
                    TXRX(idx_internal, 1), ...
                    zeros(n_TXRX_intern, 1), ...
                    TXRX(idx_internal, 2), ...
                    h_TXRX + zeros(n_TXRX_intern, 1)];

                % Remove inerior point(s) from index vector.
                idx_txrx_not_in_pd(tmp_idx(idx_in_TXRX)) = false;

                % Add all remaining points to the point_domain list.
                % Get point coordinates.
                TXRX_to_add = TXRX(idx_txrx_not_in_pd, :);
                % Create input matrix (without numbering) for above points.
                point_TXRX = [TXRX_to_add(:,1), zeros(n_TXRX_add, 1), ...
                              TXRX_to_add(:,2), h_TXRX + zeros(n_TXRX_add, 1)];
                % Add TX/RX input matrix to current point_domain list and
                % sort w.r.t. x-coordinate.
                tmp_point_domain = [sortrows([point_domain(1:end-2, 2:end); point_TXRX]); ...
                                    point_domain(end-1:end, 2:end)];
                % Prepend ongoing numbering.
                point_domain = [(1:n_domain).', tmp_point_domain];

            else
                % Set points where finer mesh is needed to be part of the
                % boundary of the model domain.

                % Insert points into domain point list.
                idx_txrx_on_domain = find(idx_txrx_not_in_pd);
                n_txrx_on_domain = length(idx_txrx_on_domain);
                n_domain = n_domain + n_txrx_on_domain;
                point_txrx_domain = [(1:n_txrx_on_domain).', ...
                                      TXRX(idx_txrx_not_in_pd, 1), ...
                                      zeros(n_txrx_on_domain, 1), ...
                                      TXRX(idx_txrx_not_in_pd, 2)];
%                                       , ...
%                                       h_TXRX + zeros(n_txrx_on_domain, 1)
                point_tmp = point_domain(end,:);
                point_domain = [point_domain(1:end-1,:); point_txrx_domain];

                % Make sure, that points will form a closed point chain
                % (w.r.t. x-coo) which forms the domain boundary outline.
                [~, idx_sort_domain] = sort(point_domain(:,2));
                point_domain = point_domain(idx_sort_domain,:);
                point_domain = [point_domain; point_tmp];
                point_domain(:,1) = 1:n_domain;
            end
        end
    end

    % Create line input ([n x 3] matrix) structure for domain boundaries.
    n_line = size(point_domain, 1);
    idx_start = 1:n_line;
    idx_end = [2:(n_line), 1];
    line = [1:n_line; idx_start; idx_end].';

    %% Create Gmsh geometry (text) file and define higher order entities.

    fileID = fopen(name, 'w');
    fprintf(fileID, '// definitions \n');
    max_off = max(TXRX(:,1))-min(TXRX(:,1));
    fprintf(fileID, 'max_off = %d;\n',max_off);
    dom_wdt = max(point_domain(:, 2))-min(point_domain(:, 2));
    fprintf(fileID, 'dom_wdt = %d;\n',dom_wdt);

    %for distance fields
    fprintf(fileID, 'lc_min = 0.05;\nlc_max = dom_wdt/2;\nDistMax = dom_wdt/0.5;\nDistMin = 0.8;\n Mesh.Algorithm=1; \n');
    %Mesh.Algorithm=1 = MeshAdapt takes a bit of time but is the least
    %likely to creates failed triangles


    % Add geometry-describing point ids.
    for i = 1:(n_domain-2)
        fprintf(fileID, 'Point(%d) = {%d, %d, %d};\n', ...
            point_domain(i,:));
    end

    for i = n_domain-1:n_domain
        fprintf(fileID, 'Point(%d) = {%d, %d, %d, dom_wdt/2};\n', ...
            point_domain(i,:));
    end


    % Add line ids.
    fprintf(fileID, '\n');
    for i = 1:n_line
        fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
            line(i,:));
    end

    % Define polygonal chain describing the complete domain boundary.
    id_line_loop = n_line + 10;
    fprintf(fileID, '\n');
    fprintf(fileID, ...
        ['Line Loop(%d) = {', repmat('%d, ', 1, n_line-1), '%d};\n'], ...
        id_line_loop, 1:n_line);

    % Define interior of this polygonal chain as closed volume (=domain).
    id_plane_surface = 10^(ceil(log10(id_line_loop)));
    fprintf(fileID, '\n');
    fprintf(fileID, ...
        'Plane Surface(%d) = {%d};\n', id_plane_surface, id_line_loop);

%     if exist('n_TXRX_intern', 'var')
%         % Add detached points.
%         fprintf(fileID, '\n');
%         for i = 1:n_TXRX_intern
%             fprintf(fileID, 'Point(%d) = {%d, %d, %d, %d};\n', ...
%                 point_internal(i,:));
%         end
%
%         % Add detached points to surface.
%         % (Only in this case they will be considered by the meshing routine.)
%         if n_TXRX_intern > 1
%             fprintf(fileID, '\n');
%             fprintf(fileID, ...
%                 ['Point{', repmat('%d, ', 1, n_TXRX_intern-1),'%d} In Surface{%d};\n'], ...
%                 point_internal(:, 1), id_plane_surface);
%
%         elseif n_TXRX_intern == 1
%             fprintf(fileID, '\n');
%             fprintf(fileID, ...
%                 'Point{%d} In Surface{%d};\n', ...
%                 point_internal(1, 1), id_plane_surface);
%         end
%     end

    % Add boundary (physical) line ids.
    % xmin (left)
    fprintf(fileID, '\n');
    fprintf(fileID, 'Physical Line("xmin") = {%d};\n', line(end, 1));
    % xmax (right)
    fprintf(fileID, 'Physical Line("xmax") = {%d};\n', line(end-2, 1));
    % ymin (surface)
    if n_line > 4
        fprintf(fileID, ['Physical Line("ymin") = {', ...
            repmat('%d, ', 1, n_line-4), '%d};\n'], ...
            line(1:(end-3), 1));
    elseif n_line == 4
        fprintf(fileID, 'Physical Line("ymin") = {%d};\n', line(1, 1));
    end
    % ymax (bottom)
    fprintf(fileID, 'Physical Line("ymax") = {%d};\n', line(end-1, 1));

    % Add physical surface id.
    if  isempty(dom_name)
        dom_name = 'res_1';
    end
    fprintf(fileID, '\n');
    fprintf(fileID, ...
        'Physical Surface("%s") = {%d};\n', dom_name, id_plane_surface);

    % Add Distance Fields
    [n_txrx,~] = size(TXRX);
    ele_idx = 2: n_txrx+1;
    for ii = 1:n_txrx
        fprintf(fileID,'Field[%d] = Distance;\n', 2*ii-1);
        fprintf(fileID,'Field[%d].NodesList = {%d};\n', 2*ii-1,ele_idx(ii));
        fprintf(fileID,'Field[%d] = Threshold;\n', 2*ii);
        fprintf(fileID,'Field[%d].IField = %d;\n', 2*ii, 2*ii-1);
        fprintf(fileID,'Field[%d].LcMin = lc_min;\n', 2*ii);
        fprintf(fileID,'Field[%d].LcMax = lc_max;\n', 2*ii);
        fprintf(fileID,'Field[%d].DistMin = DistMin;\n', 2*ii);
        fprintf(fileID,'Field[%d].DistMax = DistMax;\n \n', 2*ii);
    end

%     [n_txrx,~] = size(TXRX);
%     for i = 2:(n_txrx+1)
%         fprintf(fileID,'Field[%d] = Distance;\n', i-1);
%         fprintf(fileID,'Field[%d].NodesList = {%d};\n', (i-1),i);
%         fprintf(fileID,'Field[%d] = Threshold;\n', (i-1)*100);
%         fprintf(fileID,'Field[%d].IField = %d;\n', (i-1)*100, i-1);
%         fprintf(fileID,'Field[%d].LcMin = lc_min;\n', (i-1)*100);
%         fprintf(fileID,'Field[%d].LcMax = lc_max;\n', (i-1)*100);
%         fprintf(fileID,'Field[%d].DistMin = DistMin;\n', (i-1)*100);
%         fprintf(fileID,'Field[%d].DistMax = DistMax;\n \n', (i-1)*100);
%     end

    fprintf(fileID,'Field[%d] = Min;\n', 2*n_txrx+1);
    fprintf(fileID,'Field[%d].FieldsList = {',2*n_txrx+1);
    fprintf(fileID, [ ...
            repmat('%d, ', 1, n_txrx-1), '%d};\n'], ...
            2*(1:n_txrx));
    fprintf(fileID,'Background Field = {%d};\n', 2*n_txrx+1);


    fclose(fileID);

    if verbosity
       fprintf('done.\n');
    end

end

