function [x, w] = getQuadratureRule(order, dim, varargin)
    % Quadrature rules for reference simplex.
    %
    % Point coordinates are given in terms of two barycentric coordinates.
    % I.e. the kartesian coordinates of the 2D reference triangle /
    % simplex.
    %
    % SYNTAX
    %
    %   [x, w] = getQuadratureRule(order, dim[, varargin])
    %
    % INPUT PARAMETERS
    %
    %   order ... Scalar, order of the quadrature that should be returned.
    %   dim   ... Scalar, dim of the problem.
    %
    % OPTIONAL PARAMETERS
    %
    %   bnd  ... Vector, denoting lower and upper integration boundary (1D).
    %   type ... Char, denoting the type of 1D Gauss quadrature rule.
    %            type = {'Legendre', 'Laguerre'}
    %
    % OUTPUT PARAMETERS
    %   x ... matrix n x dim, of quadrature points.
    %   w ... vector n x 1 of corresponding weights.
    %
    % REMARKS
    %
    %   Reference simplex of dimension dim is given by vertex (cartesian) 
    %   coordinates.
    %
    %     [eye(dim), zeros(dim, 1)]
    %
    %   The triangle schemes are mostly from Strang and Fix.
    %
    % COPYRIGHT
    %   2D code originally written by Jan Blechta.
    %   1D G.-Legendre originally written by Martin Afanasjew.
    %   1D G.-Laguerre originally written by Ralph-Uwe Börner.
    
    %% Check input and set properties.
    
    % Define possible input keys and its properties checks.
    input_keys = {'bnd', 'type'};
    assertBnd = @(x) assert(isvector(x) && length(x) == 2 && x(1) < x(2), ...
        ['bnd - Vector [2 x 1] denoting the lower and the upper ', ...
        'integral boundary, expected.']);
    assertType = @(x) assert(ischar(x) && ...
        any(strcmp(x, {'Legendre', 'Laguerre'})), ...
        ['type - Char denoting quadrature type "Legendre" or ', ...
        '"Laguerre" expected.']);
    
    % Create inputParser object and set possible inputs with defaults.
    parser_obj = inputParser();
    parser_obj.addParameter(input_keys{1}, [], assertBnd);
    parser_obj.addParameter(input_keys{2}, 'Legendre', assertType);
   
    % Exctract all properties from inputParser.
    parse(parser_obj, varargin{:});
    args = parser_obj.Results;

    %% Choose propper subroutine.
    
    switch dim
        case 1
            switch args.type
                case 'Legendre'
                    [x, w] = GaussLegendre(order, args.bnd);
                case 'Laguerre'
                    [x, w] = GaussLaguerre(order);
            end
        case 2
            [x, w] = Quad2D(order);
        otherwise
            error('Only 1 and 2 dimensional quadrature rules available.');
    end

end

function [x, w] = GaussLegendre(n, bnd)
    % Provides nodes and weights of the Gauss-Legendre quadrature in 1D.
    %
    % \in_{bnd_lo}^{bnd_up} f(x) dx = ...
    %       \sum_{j = 0}^{n - 1} w_j f(x_j)
    %
    % INPUT PARAMETER
    %   n ... Positive scalar integer that denotes the number of desired
    %         nodes. The returned quadrature is exact for polynomials of
    %         degree 2 * n - 1 or less.
    %
    % OPTIONAL PARAMETER
    %   bnd ... Vector, denoting lower and upper integration boundary.
    %
    % OUTPUT PARAMETER
    %   x ... Column vector of quadrature nodes in the interval specified
    %         by 'bnd_lo' and 'bnd_up'.
    %   w ... Row vector of quadrature weights corresponding to the nodes
    %         in 'x'. The weights are rescaled to the interval specified by
    %         'bnd_lo' and 'bnd_up' such that 'sum(w) == bnd_up - bnd_lo'.
    
    % Check input
    if n < 1
        warning('Quadratur order > 0 expected. Set n = 1.');
        n = 1;
    end
    
    if isempty(bnd)
        % Set lower and upper bound w.r.t. the edge of the reference 
        % simplex.
        bnd_lo = 0;
        bnd_up = 1;
    else
        bnd_lo = bnd(1);
        bnd_up = bnd(2);
    end

    % Construct matrix for compuation of nodes and weights.
    k = (1:n - 1).';
    c = k ./ sqrt(4 .* k .^ 2 - 1);
    A = spdiags([[c; 0], [0; c]], [-1, +1], n, n);
    
    % Compute nodes and weights on [-1, 1] using eigen decomposition.
    [V, D] = eig(full(A));
    x = diag(D);
    w = V(1, :) .^ 2;
    
    % Transform nodes and weights from [-1, 1] to [bnd_lo, bnd_up].
    x = 0.5 * (bnd_lo + bnd_up) + 0.5 * (bnd_up - bnd_lo) .* x;
    w = (bnd_up - bnd_lo) .* w;
end

function [x, w] = GaussLaguerre(n)
    % Provides nodes and weights of the Gauss-Laguerre quadrature in 1D.
    %
    % \in_{0}^{\inf} \exp(-x) f(x) dx = ...
    %       \sum_{j = 0}^{n - 1} w_j f(x_j) 
    %
    % INPUT PARAMETER
    %   n ... Positive scalar integer that denotes the number of desired
    %         nodes.
    %
    % OUTPUT PARAMETER
    %   x ... Column vector of quadrature nodes in the interval 0 to \inf.
    %         The smallest abscissa is returned in x[1], the largest in 
    %         x[n].
    %   w ... Row vector of quadrature weights corresponding to the nodes
    %         in 'x'. The weights are rescaled to the interval specified by
    %         'bnd_lo' and 'bnd_up' such that 'sum(w) == bnd_up - bnd_lo'.
    %
    % For a description of the following routines see:
    % Numerical Recipes, Chapter 4, Press et al 2007

    % Set up parameters.
    tolerance = 1e-14;
    max_iter = 10;
    
    % Initiallize quantities.
    z = 0;
    x = zeros(n, 1);
    w = zeros(1, n);

    % Loop over desired roots.
    for i = 1:n
        if i == 1
            % Initial guess for the smallest root.
            z = 3 / (1 + 2.4 * n);
        elseif i == 2
            % Initial guess for the second root. 
            z = z + 15 / (1 + 2.5 * n);
        else
            % Initial guess for the other roots.
            ai = i - 2;
            z = z + (1 + 2.55 * ai) / (1.9 * ai) * (z - x(ai));
        end
        
        % Refinement by Newton's method.
        for its = 1:max_iter
            p1 = 1;
            p2 = 0;
            
            % Loop up the recurrence relation to get the Laguerre 
            % polynomial evaluated at z.
            for j = 1:n
                p3 = p2;
                p2 = p1;
                p1 = ((2 * j - 1 - z) * p2 - (j - 1) * p3) / j;
            end
            
            % p1 is now the desired Laguerre polynomial. 
            % Now compute pp, its derivative, by a standard relation 
            % involving also p2, the polynomial of one lower order.
            pp = n * (p1 - p2) / z;
            z1 = z;
            
            % Apply Newton’s formula.
            z  = z1 - p1 / pp;
            if(abs(z - z1) <= tolerance) 
                break;
            end
        end
        x(i) = z;
        w(i) = -1 / (pp * n * p2);
    end
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
        case 3
            % Scheme from Strang and Fix, 6 points, degree of precision 3
            data = [
                0.659027622374092, 0.231933368553031, 1.0/12.0;
                0.659027622374092, 0.109039009072877, 1.0/12.0;
                0.231933368553031, 0.659027622374092, 1.0/12.0;
                0.231933368553031, 0.109039009072877, 1.0/12.0;
                0.109039009072877, 0.659027622374092, 1.0/12.0;
                0.109039009072877, 0.231933368553031, 1.0/12.0;
            ];
        case 4
            % Scheme from Strang and Fix, 6 points, degree of precision 4
            data = [
                0.816847572980459, 0.091576213509771, 0.109951743655322/2.0;
                0.091576213509771, 0.816847572980459, 0.109951743655322/2.0;
                0.091576213509771, 0.091576213509771, 0.109951743655322/2.0;
                0.108103018168070, 0.445948490915965, 0.223381589678011/2.0;
                0.445948490915965, 0.108103018168070, 0.223381589678011/2.0;
                0.445948490915965, 0.445948490915965, 0.223381589678011/2.0;
            ];
        case 5
            % Scheme from Strang and Fix, 7 points, degree of precision 5
            data = [
                0.33333333333333333, 0.33333333333333333, 0.22500000000000000/2.0;
                0.79742698535308720, 0.10128650732345633, 0.12593918054482717/2.0;
                0.10128650732345633, 0.79742698535308720, 0.12593918054482717/2.0;
                0.10128650732345633, 0.10128650732345633, 0.12593918054482717/2.0;
                0.05971587178976981, 0.47014206410511505, 0.13239415278850616/2.0;
                0.47014206410511505, 0.05971587178976981, 0.13239415278850616/2.0;
                0.47014206410511505, 0.47014206410511505, 0.13239415278850616/2.0;
            ];
        case 6
            % Scheme from Strang and Fix, 12 points, degree of precision 6
            data = [
                0.873821971016996, 0.063089014491502, 0.050844906370207/2.0;
                0.063089014491502, 0.873821971016996, 0.050844906370207/2.0;
                0.063089014491502, 0.063089014491502, 0.050844906370207/2.0;
                0.501426509658179, 0.249286745170910, 0.116786275726379/2.0;
                0.249286745170910, 0.501426509658179, 0.116786275726379/2.0;
                0.249286745170910, 0.249286745170910, 0.116786275726379/2.0;
                0.636502499121399, 0.310352451033785, 0.082851075618374/2.0;
                0.636502499121399, 0.053145049844816, 0.082851075618374/2.0;
                0.310352451033785, 0.636502499121399, 0.082851075618374/2.0;
                0.310352451033785, 0.053145049844816, 0.082851075618374/2.0;
                0.053145049844816, 0.636502499121399, 0.082851075618374/2.0;
                0.053145049844816, 0.310352451033785, 0.082851075618374/2.0;
            ];
        otherwise
            error('Order %d unsupported.', order);
    end
    x = data(:, 1:2);
    w = data(:, 3);
end