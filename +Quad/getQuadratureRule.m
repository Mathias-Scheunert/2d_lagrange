function [x, w] = getQuadratureRule(order, dim)
    % Quadrature rules for reference simplex.
    %
    % Point coordinates are given in terms of two barycentric coordinates.
    % I.e. the kartesian coordinates of the 2D reference triangle /
    % simplex.
    %
    % SYNTAX
    %
    %   [x, w] = getQuadratureRule(order, dim)
    %
    % INPUT PARAMETERS
    %
    %   order ... Scalar, order of the quadrature that should be returned.
    %   dim   ... Scalar, dim of the problem.
    %
    % OUTPUT PARAMETERS
    %   x ... matrix n x dim, of quadrature points.
    %   w ... vector n x 1 of corresponding weights.
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
    %   2D code originally written by Jan Blechta (CurlCurl-Toolbox).
    %   1D code originally written by Martin Afanasjew (toolbox).
    
    switch dim
        case 1
            [x, w] = Quad1D(order);
        case 2
            [x, w] = Quad2D(order);
        otherwise
            error('Only 1 and 2 dimensional quadrature rules available.');
    end

end

function [x, w] = Quad1D(n)
    % Provides nodes and weights of the Gauss-Legendre quadrature in 1D.
    %
    % OUTPUT PARAMETER
    %   x ... Column vector of quadrature nodes in the interval specified
    %         by 'bnd_lo' and 'bnd_up'.
    %   w ... Row vector of quadrature weights corresponding to the nodes
    %         in 'x'. The weights are rescaled to the interval specified by
    %         'bnd_lo' and 'bnd_up' such that 'sum(w) == bnd_up - bnd_lo'.
    %
    % INPUT PARAMETER
    %   n ... Positive scalar integer that denotes the number of desired
    %         nodes. The returned quadrature is exact for polynomials of
    %         degree 2 * n - 1 or less.
    %         % TODO: As the code is currently designed for 1D and 2D only
    %         this number corrensponds to the dim of the problem.
    
    % Set lower and upper bound w.r.t. the edge of the reference simplex.
    bnd_lo = 0;
    bnd_up = 1;

    % Construct matrix for compuation of nodes and weights.
    k = (1:n - 1).';
    c = k ./ sqrt(4 .* k .^ 2 - 1);
    A = spdiags([[c; 0], [0; c]], [-1, +1], n, n);
    
    % Compute nodes and weights on [-1, 1] using eigen decomposition.
    [V, D] = eig(full(A));
    x = diag(D);
    w = V(1, :) .^ 2;
    
    % Transform nodes and weights from [-1, 1] to [a, b].
    x = 0.5 * (bnd_lo + bnd_up) + 0.5 * (bnd_up - bnd_lo) .* x;
    w = (bnd_up - bnd_lo) .* w;
end

function [x, w] = Quad2D(order)
    % Provides nodes and weights of the Gauss-Legendre quadrature in 2D.
    %
    % OUTPUT PARAMETER
    %   x ... Column vector of quadrature nodes in the interval specified
    %         by 'bnd_lo' and 'bnd_up'.
    %   w ... Row vector of quadrature weights corresponding to the nodes
    %         in 'x'. The weights are rescaled to the interval specified by
    %         'bnd_lo' and 'bnd_up' such that 'sum(w) == bnd_up - bnd_lo'.
    %
    % INPUT PARAMETER
    %   oder ... Positive scalar, denoting the order of the Lagrange
    %            elements.
    
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