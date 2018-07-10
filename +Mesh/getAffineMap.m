function map = getAffineMap(cell_num, mesh, point)
    % Provides barycentric coordinate transformation for a point in a cell.
    %
    % The linear map matrix and translation vector for the barycentric
    % coordinate system representation within a cell cell_num of a point, 
    % given in the cartesian coordinates, is calculated.
    % 
    % SYNTAX
    %   map = getAffineMap(cell_num, mesh[, point])
    %
    % INPUT PARAMETER
    %   cell_num ... Scalar, denoting the cell number.
    %   mesh     ... Struct, containing mesh information, i.e. coordinates
    %                of vertices and its relation to the triangles and 
    %                edges.
    %
    % OPTIONAL PARAMETER
    %   point ... Matrix nx2, containing coordinates of observation points.
    %
    % OUTPUT PARAMETER
    %   map ... Struct, containing linear map matrix (and it's inverse) and
    %           the translation vector required to map between the
    %           baryzentric coordinate system representation and the
    %           cartesian coordinate system representation (and back) 
    %           w.r.t. a specific cell (and of a given point).
    %
    % REMARKS
    %
    %   Consider a Point p = [x; y] within a triangle with the three 
    %   vertices r_1 = [x_1; y_1], r_2 = [x_2; y_2], r_3 = [x_3; y_3]
    %
    %   For the barycentric coordinates lambda_1, lambda_2, lambda_3 it
    %   holds:
    %   [lambda_1; lambda_2] = T^(-1) * (p - r_3)
    %   lambda_3             = 1 - lambda_1 - lambda_2
    %   and
    %   p = T * [lambda_1; lambda_2] + r_3
    %   with
    %   T = [(x_1 - x_3) , (x_2 - x_3);
    %        (y_1 - y_3) , (y_2 - y_3)]
    %   
    %   To avoid forming the inverse (even using \ is quite slow) by
    %   lambda_1 = (y_2 - y_3)*(x - x_3) + (x_3 - x_2)*(y - y_3) / det(T)
    %   lambda_2 = (y_3 - y_1)*(x - x_3) + (x_1 - x_3)*(y - y_3) / det(T)
    %   lambda_3 = 1 - lambda_1 - lambda_2
    %   
    %   Alternatively one can extent the expressions as given by
    %   (explicitly include the redundand lambda_3)
    %   [lambda_1; lambda_2; lambda_3] = R^(-1) * [p; 1]
    %   and 
    %   [p; 1] = R * [lambda_1; lambda_2; lambda_3]
    %   with
    %   R = [x_1, x_2, x_3;
    %        y_1, y_2, y_3;
    %        1  , 1  , 2   ]
    %
    %   Resulting coordinates in reference simplex:
    %   -> lambda_1 = x_hat, lambda_2 = y_hat
    %
    %   local
    %
    %   y_hat
    %   ^
    %   '-> x_hat
    % 
    %   [0,1]=point 3
    %   |    \
    %   |     \
    %   3^     2v
    %   |       \
    %   |        \
    %    - - 1> - - [1,0]=point 2  
    %   [0,0]=point 1
    %
    %   With the definitions for B and b: 
    %      reference             global
    %       point 1   maps to  [x_3,y_3]
    %       point 2   maps to  [x_1,y_1]
    %       point 3   maps to  [x_2,y_2]
    %   And definitions from Mesh.appendElementInfo.m
    %       edge 1    maps to   edge 3
    %       edge 2    maps to   edge 1
    %       edge 3    maps to   edge 2

    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    assert(isscalar(cell_num) && cell_num <= length(mesh.cell2vtx), ...
        'cell_num - scalar <= numbers of cells expected.');
    if nargin == 3
        assert(ismatrix(point) && size(point, 2) == 2, ...
            'point - matrix nx2, containing coordinates, expected.');
    end
    
    %% Prepare mapping.
    
    % Get coordinates of cell vertices.
    vert_cord = mesh.cell2cord{cell_num};
    
    %% Create mapping from barycentric to cartesian coordinates.
    
    % Define linear mapping.
    % (neglecting redundand third barycentric coordinate)
    B = [vert_cord(1,1) - vert_cord(3,1), vert_cord(2,1) - vert_cord(3,1); ...
         vert_cord(1,2) - vert_cord(3,2), vert_cord(2,2) - vert_cord(3,2)];

%     % (explicitly incorporate third barycentric coordinate)
%     B_long = [vert_cord(1,1), vert_cord(2,1), vert_cord(3,1); ...
%          vert_cord(1,2), vert_cord(2,2), vert_cord(3,2); ...
%          1             , 1             , 1];
     
    % Define translation part.
    % (neglecting redundand third barycentric coordinate)
    b = [vert_cord(3,1); vert_cord(3,2)];
    
%     % (explicitly incorporate third barycentric coordinate)
%     b_long = [vert_cord(3,1); vert_cord(3,2); 1];

    % Define transposed invers by using the adjoint matrix.
    BinvT = (1 / det(B)) * ([1, -1; -1, 1] .* rot90(B, 2));
    
    % Set association between the index of the
    % local/reference simplex vertices to the global/mapped ones.
    loc2glo = [3, 1, 2];
    
    % Summarize infos.
    map = struct();
    map.B = B;
    map.b = b;
    map.detB = det(B);
    map.BinvT = BinvT;
    map.loc2glo = loc2glo;
    
    %% Create mapping from cartesian to barycentric coordinates.
 
    if nargin == 3
        
        % Use short explicit formulation from wikipadia.org. to obtain the
        % two barycentric coordinates (coordinates in reference simplex)
        % w.r.t. the given point(s).
        map = @(input) [(vert_cord(2,2) - vert_cord(3,2)) * (input(1) - vert_cord(3,1)) ...
                + (vert_cord(3,1) - vert_cord(2,1)) * (input(2) - vert_cord(3,2)); ...
                  (vert_cord(3,2) - vert_cord(1,2)) * (input(1) - vert_cord(3,1)) ...
                + (vert_cord(1,1) - vert_cord(3,1)) * (input(2) - vert_cord(3,2))] ...
                * 1/map.detB; 
        if size(point, 1) == 1
            lambda = map(point).';
        else
            point = mat2cell(point, ones(size(point, 1), 1), 2);
            lambda = cell2mat(cellfun(@(x) {map(x).'}, point));
        end
        
        % Summarize infos.
        map = struct();
        map.xy_ref = lambda;
    end
end