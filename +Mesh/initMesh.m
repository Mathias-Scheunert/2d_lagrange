function mesh = initMesh(var, varargin)
    % Creates and optionaly uniformly refines a 2D triangle-based mesh.
    %
    % SYNTAX
    %   mesh = initMesh(var, varargin)
    %
    % INPUT PARAMETER
    %   var ... Char, denoting the variant of elementary grid.
    %
    % OPTIONAL PARAMETER
    %   bnd ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   TX        ... Vector [n x 2], denoting the source position(s).
    %   RX        ... Vector [m x 2], denoting the receiver position(s).
    %   topo      ... Vector [k x 2], denoting the descrete topography.
    %   ref       ... Scalar, denoting the number of uniform refinement
    %                 steps.
    %   verbosity ... Logical, denoting if current status should be
    %                 printed
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges
    %            as well as the relation between triangles and edges.
    %
    % Elemental information contained in mesh:
    % (Note that further relations can be obtained by append... routines in
    %  +Mesh folder)
    % dim      ... Scalar, denoting the geometrical dimension.
    % vertices ... Matrix [n x 2], containing the x,y coordinates the n
    %              vertices in mesh.
    % cell2vtx ... Matrix [m x 3], containting the list of vertices each of
    %              the m cells is formed of, respectively.
    % parameter_domain      ... Vector [m x 1], containing k different
    %                           parameter domain indices, each of the m
    %                           cells is related to.
    % parameter_domain_name ... Cell of chars [k x 1], containing the
    %                           names/tags for each of the k parameter
    %                           domains.
    % bnd_edge_part         ... Vector [l x 1], containing p different
    %                           boundary indices, each of the l boundary
    %                           edges is related to.
    % bnd_edge_part_name    ... Cell of chars [p x 1], containing the
    %                           names/tags for each of the p boundary
    %                           parts.
    % bnd_edge2vtx          ... Matrix [l x 2], containting the list of
    %                           vertices each of l boundary edges is
    %                           formed of, respectively.
    % -> To be able to use external meshes for solving PDEs, these
    %    information sould be provided.

    %% Check input.

    assert(ischar(var), ...
        'var - Char, denoting the basic grid type, expected.');

    % Define possible input keys and its properties checks.
    input_keys = {'bnd', 'ref', 'TX', 'RX', 'topo', 'dom_name', 'name', ...
                  'dike', 'verbosity'};
    assertBnd = @(x) assert(isvector(x) && length(x) == 4 && ...
        x(1) < x(2) && x(3) < x(4), ...
        ['bnd - Vector [1 x 4] denoting the lower and the upper ', ...
        'integral boundary, expected.']);
    assertRef = @(x) assert(isscalar(x) && ~islogical(x) && x >= 0, ...
        'ref - Scalar, denoting the number of uniform ref steps, expected.');
    assertPos = @(x) assert(isempty(x) || (ismatrix(x) && size(x, 2) == 2), ...
        'TX/RX/topo - Vector [n x 2], denoting positions, expected.');
    assertName = @(x) assert(ischar(x), ...
        'name - Character of Gmsh file name to load (including file extention!).');
    assertDomName = @(x) assert(ischar(x), ...
        'dom_name - Character of parameter domain name, expected.');
    assertDike = @(x) assert(isvector(x) && length(x) == 2, ...
        'dike - Vector [1 x 2] of left boundary and width,m expected.');
    assertVerbose = @(x) assert(islogical(x), ...
        'verbosity - logical, denoting if status should be printed, expected');

    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, [-1, 1, -1, 1], assertBnd);
    parser_obj.addParameter(input_keys{2}, 0, assertRef);
    parser_obj.addParameter(input_keys{3}, [], assertPos);
    parser_obj.addParameter(input_keys{4}, [], assertPos);
    parser_obj.addParameter(input_keys{5}, [], assertPos);
    parser_obj.addParameter(input_keys{6}, [], assertName);
    parser_obj.addParameter(input_keys{7}, 'def', assertDomName);
    parser_obj.addParameter(input_keys{8}, [5, 5], assertDike);
    parser_obj.addParameter(input_keys{9}, false, assertVerbose);

    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    %% Create mesh.

    switch var
        case 'cube'
            mesh = Mesh.createUnitCubeMesh(args.bnd, [1, 1], args.verbosity);
        case 'disc'
            mesh = Mesh.createUnitDiscMesh(args.ref, args.verbosity);
            args.ref = 0;
        case 'gmsh_dike'
            mesh = Mesh.createVerticalDikeMesh(args.bnd, args.TX, ...
                                               args.dike(1), args.dike(2), ...
                                               args.ref);
        case 'gmsh_create'
            mesh = Mesh.createGmsh(args.bnd, args);
        case 'gmsh_load'
            if isempty(args.name)
                error('Provide file name as an argument for Mesh.initMesh().');
            end
            mesh = Mesh.loadGmsh(args.name, ...
                'verbosity', args.verbosity, ...
                'ref', args.ref);
        otherwise
            error('Unknown mesh type.');
    end

    % Append edge information.
    if args.verbosity
       fprintf('Append edge info ... ');
    end
    mesh = Mesh.appendEdgeInfo(mesh);
    if args.verbosity
       fprintf('done.\n');
    end

    % Append mapping between information and coordinates.
    if args.verbosity
       fprintf('Append Coo map ... ');
    end
    mesh = Mesh.appendCoordInfo(mesh);
    if args.verbosity
       fprintf('done.\n');
    end

    % Append boundary information.
    if args.verbosity
       fprintf('Append BND info ... ');
    end
    mesh = Mesh.appendBndInfo(mesh);
    if args.verbosity
       fprintf('done.\n');
    end

    %% Refine mesh.

    if any(strcmp(mesh.type, {'cube'}))
        if args.verbosity
            fprintf('Refine mesh ... ');
        end
        mesh = Mesh.refineMeshUniform(mesh, args.ref);
        if args.verbosity
           fprintf('done.\n');
        end
    end

    if args.verbosity
        fprintf('... Mesh struct initialized.\n \n');
    end

    %% Get affine mappings.

    if args.verbosity
       fprintf('Obtain affine maps ... ');
    end
    % Get maps.
    mesh.maps = cellfun(@(x) {Mesh.getAffineMap(x, mesh)}, ...
                    num2cell(1:size(mesh.cell2vtx, 1)).');

    % Exclude vertex and edge relations.
    mesh.loc2glo = mesh.maps{1}.loc2glo;
    [~, mesh.glo2loc] = ismember(1:3, mesh.loc2glo);
    mesh.loc_normals = mesh.maps{1}.loc_normals;
    mesh.loc2glo_orient = mesh.maps{1}.loc2glo_edge_dir;
    if args.verbosity
       fprintf('done.\n');
    end
end
