function mesh = createUnitCubeMesh(bnd, divisions, verbosity)
    % Build unit cube mesh.
    %
    % SYNTAX
    %   mesh = createUnitCubeMesh(bnd, divisions[, verbosity])
    %
    % INPUT PARAMETER
    %
    %   bnd       ... Boundaries of modeling area [xmin, xmax, ymin, ymax].
    %   divisions ... Number of divisions along axes, i.e., [nx, ny].
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if verbose output is desired.
    %
    % OUTPUT PARAMETER
    %   mesh      ... Structure with basic topoogy and geometry data.
    %
    % REMARKS
    %
    %   generate_unit_cube_mesh([3, 5]) results in unit square mesh of
    %   triangles with 3 divisiones along x and 5 along y axis. Analogically,
    %   generate_unit_cube_mesh([3, 4, 5]) gives mesh of tetrahedrons.
    %
    %   The mesh is a struct with arrays vertex_coords and cells.
    %   Value of vertex_coords(i, j) is i-th coordinate of j-th vertex.
    %   Value cells(:, j) gives vertex indices constituting j-th cell.
    %
    % COPYRIGHT
    %   Code originally written by Jan Blechta (CurlCurl-Toolbox).
    
    %% Check input.
    
    assert(isvector(bnd) && length(bnd) == 4, ...
        'bnd - Vector [4 x 1] denoting the domain boundaries, expected.');
    assert(isvector(divisions) && length(divisions) == 2, ...
        ['divisions - Vector [2 x 1] denoting the number of divisions ', ...
        'along axes, expected.']);
    
    if nargin < 3
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, denoting if verbose output is desired, expected.');        
    end
    
    %% Create mesh.
    
    if verbosity
       fprintf('Create basic mesh ... \n'); 
    end

    mesh = generate_unit_2cube_mesh(bnd, divisions(1), divisions(2));
    
    if verbosity
       fprintf('done.\n');
    end
end

function mesh = generate_unit_2cube_mesh(bnd, nx, ny)
    % Ensure consistency.
    if nx == 0
        nx = 1;
    end
    if ny == 0 
        ny = 1;
    end

    % Build vertex coordinates.
    x = linspace(bnd(1), bnd(2), nx + 1);
    y = linspace(bnd(3), bnd(4), ny + 1);
    [x, y] = ndgrid(x, y);
    num_vertices = (nx + 1)*(ny + 1);
    vertex_coords = zeros(2, num_vertices);
    vertex_coords(1, :) = reshape(x, num_vertices, 1);
    vertex_coords(2, :) = reshape(y, num_vertices, 1);

    % Build cell-to-vertex topology.
    cells = zeros(3, 2*nx*ny);
    for iy = 0:ny-1
        for ix = 0:nx-1
            v0 = 1 + iy*(nx + 1) + ix;
            v1 = v0 + 1;
            v2 = v0 + (nx + 1);
            v3 = v1 + (nx + 1);

            c0 = 1 + 2*(iy*nx + ix);
            cells(:, c0+0) = [v0, v1, v2];
            cells(:, c0+1) = [v1, v2, v3];
        end
    end

    % Summarize.
    mesh = struct();
    mesh.dim = 2;    
    mesh.vertices = vertex_coords.';
    mesh.cell2vtx = cells.';
    mesh.parameter_domain = ones(size(mesh.cell2vtx, 1), 1);
    mesh.parameter_domain_name = {'entire'};
end