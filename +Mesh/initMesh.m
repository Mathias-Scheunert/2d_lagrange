function mesh = initMesh(var, bnd, ref, verbosity)
    % Creates and optionaly uniformly refines a 2D triangle-based mesh.
    % 
    % SYNTAX
    %   mesh = initMesh(bnd, var, ref[, verbosity])
    %
    % INPUT PARAMETER
    %   var ... Char, denoting the variant of elementary grid.
    %   bnd ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   ref ... Scalar, denoting the number of uniform refinement steps.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %               of vertices and its relation to the triangles and 
    %               edges as well as the relation between triangles and 
    %               edges.
    
    %% Check input.
    
    assert(ischar(var), ...
        'var - Char, denoting the basic grid type, expected.');
    assert(isscalar(ref) && ref >= 0, ...
        'ref - Scalar, denoting the number of uniform ref steps, expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    %% Create mesh.
    
    if verbosity
       fprintf('Create basic mesh ... '); 
    end
    switch var
        case 'rhomb'
            mesh = Mesh.createRhombMesh(bnd);
        case 'cube'
            mesh = Mesh.createUnitCubeMesh(bnd, [3, 3]);
        otherwise 
            error('Unknown mesh type.');
    end
    if verbosity
       fprintf('done.\n'); 
    end
    mesh.type = var;
    
    % Append edge information.
    if verbosity
       fprintf('Append edge info ... '); 
    end
    mesh = Mesh.appendElementInfo(mesh);
    if verbosity
       fprintf('done.\n'); 
    end
    
    % Append mapping between information and coordinates.
    if verbosity
       fprintf('Append Coo map ... '); 
    end
    mesh = Mesh.appendMappings(mesh);
    if verbosity
       fprintf('done.\n'); 
    end
    
    % Append boundary information.
    if verbosity
       fprintf('Append BND info ... '); 
    end
    mesh.bnd = bnd;
    mesh = Mesh.appendBndInfo(mesh);
    if verbosity
       fprintf('done.\n'); 
    end
    
    %% Refine mesh.

    if verbosity
        fprintf('Refine mesh ... '); 
    end
    mesh = Mesh.refineMeshUniform(mesh, ref);
    if verbosity
       fprintf('done.\n');
       fprintf('... Mesh struct initialized.\n \n'); 
    end
    
    %% Get affine mappings.
    
    if verbosity
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
    if verbosity
       fprintf('done.\n'); 
    end
end

