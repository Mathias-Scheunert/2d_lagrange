function mesh = createRhombMesh(bnd, verbosity)
    % Creates a mesh in 2D representing a rhombus in a finit rectangle.
    % 
    % SYNTAX
    %   mesh = createRhombMesh(bnd[, verbosity])
    %
    % INPUT PARAMETER
    %   bnd ... Vector, denoting boundaries of modeling area  
    %           [xmin, xmax, ymin, ymax].
    %
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing vertex2coordinates (vertices), 
    %            simplex2vertex (cell2vtx), and simplex parameter domains 
    %            (parameter_domains).
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if verbose output is desired.
    %
    % REMARKS
    %
    %   The points are referred to a coordinate system using x(1), y(1) as
    %   coordinate origin.
    %   The lexicographic order is from left to right and from bottom to 
    %   top.

    %% Check input.
    
    assert(isvector(bnd) && length(bnd) == 4, ...
        'Expected bnd to be vectors 4 x 1 denoting the domain boundaries.');
    
    if nargin < 2
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - Logical, denoting if verbose output is desired, expected.');        
    end

    %% Create vertex (point) list.
    
    if verbosity
       fprintf('Create basic mesh ... '); 
    end
    
    % Line number = number of vertex
    % 1. colum = x coordinate
    % 2. colum = y coordinate
    
    % Construct regular grid as basis for unstructured grid.
    n = 5;
    x_bnd_node = linspace(bnd(1), bnd(2), n);
    y_bnd_node = linspace(bnd(3), bnd(4), n);
    
    % Collect boundary vertices.
    x_bnd = x_bnd_node.';
    y_bnd = y_bnd_node(y_bnd_node ~= mean(y_bnd_node)).';
    vert_bnd = [[x_bnd, bnd(3) + zeros(n, 1)];
                [x_bnd, bnd(4) + zeros(n, 1)];
                [bnd(1) + zeros(n-3, 1), y_bnd(2:end-1)];
                [bnd(2) + zeros(n-3, 1), y_bnd(2:end-1)]];
            
    % Collect vertices from rhombus.
    vert_rhomb = [x_bnd_node((n + 1)/2), y_bnd_node((n + 3)/4);
                  x_bnd_node([(n + 3)/4, (3*n + 1)/4]).', ...
                    y_bnd_node((n + 1)/2) + zeros(2, 1); ...
                  x_bnd_node((n + 1)/2), y_bnd_node((3*n + 1)/4)];
  
    % Summarize vertices and sort them with respect to local coordinate
    % system in lexicographic order.
    vert_list = sortrows([vert_bnd; vert_rhomb], [2, 1]);

    %% Create cell (triangle) list.
    
    % Line number = number of cell.
    % colums = 1st, 2nd, 3rd vertex forming a cell.
    
    % Sorted by a lexocigraphic order of the triangles center.
    cell_list = [1, 2, 6;
                 2, 3, 7;
                 2, 6, 9;
                 2, 7, 9;
                 3, 4, 7;
                 4, 5, 8;
                 4, 7, 10;
                 4, 8, 10;               
                 6, 9, 11;
                 7, 9, 12;
                 7, 10, 12;
                 8, 10, 13;
                 9, 11, 15;
                 9, 12, 15;
                 10, 12, 17;
                 10, 13, 17;
                 11, 14, 15;
                 12, 15, 16;
                 12, 16, 17;
                 13, 17, 18];
                 
    %% Summarize information.
    
    mesh = struct();
    mesh.dim = 2;
    mesh.vertices = vert_list;
    mesh.cell2vtx = cell_list;
    mesh.parameter_domain = ones(size(mesh.cell2vtx, 1), 1);
    mesh.parameter_domain_name = {'entire'};
    
    if verbosity
       fprintf('done.\n');
    end
end
