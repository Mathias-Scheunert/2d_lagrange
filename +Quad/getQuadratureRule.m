function [x, w] = getQuadratureRule(order)
    % Quadrature rules for reference simplex.
    %
    % Point coordinates are given in terms of two barycentric coordinates.
    % I.e. the kartesian coordinates of the 2D reference triangle /
    % simplex.
    %
    % SYNTAX
    %
    %   [x, w] = getQuadratureRule(order)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   order ... Order of the quadrature that should be returned.
    %   x     ... An n x dim matrix of quadrature points.
    %   w     ... An n x 1 vector of corresponding weights.
    %
    % REMARKS
    %
    %   Reference simplex of dimension dim is given by
    %   vertex (cartesian) coordinates.
    %
    %     [eye(dim), zeros(dim, 1)]
    %
    %   The triangle schemes are mostly from Strang and Fix.
    %
    % COPYRIGHT
    %   Code originally written by Jan Blechta (CurlCurl-Toolbox).
    
    switch order
        case {0, 1}
            % Scheme from Zienkiewicz and Taylor, 1 point, degree of precision 1
            data = [
                1.0/3.0, 1.0/3.0, 1.0/2.0;
            ];
        case 2
            % Scheme from Strang and Fix, 3 points, degree of precision 2
            data = [
                1.0/6.0, 1.0/6.0, 1.0/6.0;
                1.0/6.0, 2.0/3.0, 1.0/6.0;
                2.0/3.0, 1.0/6.0, 1.0/6.0;
            ];
        otherwise
            error('Order %d unsupported.', order);
    end
    
    x = data(:, 1:2);
    w = data(:, 3);
end