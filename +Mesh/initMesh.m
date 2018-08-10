function mesh = initMesh(var, bnd, varargin)
    % Creates and optionaly uniformly refines a 2D triangle-based mesh.
    % 
    % SYNTAX
    %   mesh = initMesh(bnd, var, varargin)
    %
    % INPUT PARAMETER
    %   var ... Char, denoting the variant of elementary grid.
    %   bnd ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %
    % OPTIONAL PARAMETER
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
    
    %% Check input.
    
    assert(ischar(var), ...
        'var - Char, denoting the basic grid type, expected.');
    assert(isvector(bnd) && length(bnd) == 4 && ...
        bnd(1) < bnd(2) && bnd(3) < bnd(4), ...
        ['bnd - Vector [1 x 4] denoting the lower and the upper ', ...
        'integral boundary, expected.']);
    
    % Define possible input keys and its properties checks.
    input_keys = {'ref', 'TX', 'RX', 'topo', 'verbosity'};
    assertRef = @(x) assert(isscalar(x) && ~islogical(x) && x >= 0, ...
        'ref - Scalar, denoting the number of uniform ref steps, expected.');
    assertPos = @(x) assert(ismatrix(x) && size(x, 2) == 2, ...
        'TX/RX/topo - Vector [n x 2], denoting positions, expected.');
    assertVerbose = @(x) assert(islogical(x), ...
        'verbosity - logical, denoting if status should be printed, expected');
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, 0, assertRef);
    parser_obj.addParameter(input_keys{2}, [], assertPos);
    parser_obj.addParameter(input_keys{3}, [], assertPos);
    parser_obj.addParameter(input_keys{4}, [], assertPos);
    parser_obj.addParameter(input_keys{5}, false, assertVerbose);
   
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;
    
    %% Create mesh.
    
    if args.verbosity
       fprintf('Create basic mesh ... '); 
    end
    switch var
        case 'rhomb'
            mesh = Mesh.createRhombMesh(bnd);
        case 'cube'
            mesh = Mesh.createUnitCubeMesh(bnd, [3, 3]);
        case 'gmsh_create'
            mesh = Mesh.createGmsh(bnd, args);
        case 'gmsh_load'
            error('Not implemented yet.');
            % TODO: somehow exclude (without poluting the repo) getGmshOutput
        otherwise 
            error('Unknown mesh type.');
    end
    if args.verbosity && ~any(strcmp(var, {'gmsh_create', 'gmsh_load'}))
       fprintf('done.\n'); 
    end
    mesh.type = var;
    
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
    mesh.bnd = bnd;
    mesh = Mesh.appendBndInfo(mesh);
    if args.verbosity
       fprintf('done.\n'); 
    end
    
    %% Refine mesh.

    if any(strcmp(mesh.type, {'cube', 'rhomb'}))
        if args.verbosity
            fprintf('Refine mesh ... '); 
        end
        mesh = Mesh.refineMeshUniform(mesh, args.ref);
        if args.verbosity
           fprintf('done.\n');
           fprintf('... Mesh struct initialized.\n \n'); 
        end
    end
    
    %% Get affine mappings.
    
    if args.verbosity
       fprintf('Obtain affine maps ... '); 
    end
    % Get maps.
    mesh.maps = cellfun(@(x) {Mesh.getAffineMap(x, mesh)}, ...
                    num2cell(1:length(mesh.cell2vtx)).');
                
    % Check consistency.
	maps = [mesh.maps{:}].';
    maps = reshape([maps(:).loc2glo], ...
        size(maps(1).loc2glo, 2), size(maps, 1)).';
    all((maps(:,1) == maps(:,:)));
    assert(size(unique(maps, 'rows'), 1) == 1, ...
        'Affine maps have different local to global vertex relations.');
    
    % Exclude vertex relations.
    mesh.loc2glo = maps(1,:);
    if args.verbosity
       fprintf('done.\n'); 
    end
end