function ele2L_fun = getElectrodeAt2L(rho1, rho2, h, I, TX)
    % Get semi-analytic solution for point source at 2 layered half-space.
    %
    % Provides relevant information of the potential of a point source
    % above a 2-layered half-sapce in 3D:
    %   \Phi(r) = (\rho1 I) / (2 \pi \abs(r - r')) ...
    %                * [1 + 2 * \sum_{j=1}^{\infty} R12^j ...
    %                / (\sqrt(1 + 4 * (j^2 h^2 / \abs(r - r')^2))]
    %         r = observation point = \sqrt(x^2 + y^2 + z^2)
    %         r'= source position
    %         R12 = (rho2 - rho1) / (rho1 + rho2)
    %
    % SYNTAX
    %   ele2L_fun = getElectrodeAt2L(rho1, rho2, h, I, TX)
    %
    % INPUT PARAMETER
    %   rho12 ... Scalar, denoting the layer resistivities.
    %   h     ... Scalar, denoting the first layer thickness.
    %   I     ... Scalar, denoting the source current.
    %   TX    ... Vector [2, x 1], denoting the source position.
    %
    % OUTPUT PARAMETER
    %   ele2L_fun ... Struct, containing function handle.

    %% Check input.

    assert(isscalar(rho1), 'rho1 - Scalar denoting resistivity, expected.');
    assert(isscalar(rho2), 'rho2 - Scalar denoting resistivity, expected.');
    assert(isscalar(h), 'h - Scalar denoting first layer thickness, expected.');
    assert(isscalar(I), 'I - Scalar denoting source current, expected.');
    assert(isvector(TX) && length(TX) == 2, ...
        'TX - Vector [2 x 1], denoting source position, expected.');

    %% Define function.

    R12 = (rho2 - rho1) / (rho1 + rho2);
    r = @(X, Y) norm([X; Y] - TX(:));
    ele2L_fun.f = @(X, Y) (rho1 * I) / (2 * pi * r(X, Y)) * ...
                          (1 + 2 * get_summation(r(X, Y), h, R12));

    %% Define gradient.

    % Skip derivation.

    %% Set required quadrature order.

    ele2L_fun.quad_ord = 4;
end

function sum = get_summation(r, h, R12)
    % Approximate the infinite summation.

    % Set up problem.
    max_iter = 1000;
    sum = R12 / sqrt(1 + (4 * h^2 / r.^2)); % initial summand
    tol = 1e-5 * norm(sum);

    % Loop until convergence is reached.
    for ii = 2:max_iter
       sum_add = R12^ii / sqrt(1 + (4 * ii^2 * h^2 / r.^2));
       if norm(sum_add) < tol
           break;
       end
       sum = sum + sum_add;
    end
end
