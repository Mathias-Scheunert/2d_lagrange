function eleVS_fun = getElectrodeAtVS(rho, I, TX)
    % Get information for an point source in a homogeneous full-space.
    %
    % Provides relevant information of the potential of a point source
    % in a homogeneous full-space in 3D:
    %   \Phi(r) = (\rho I) / (4 \pi \abs(r - r'))
    %         r = observation point = \sqrt(x^2 + y^2 + z^2)
    %         r'= source position
    %
    % SYNTAX
    %   eleVS_fun = getElectrodeAtVS(rho, I, TX)
    %
    % INPUT PARAMETER
    %   rho ... Scalar, denoting the VS resistivity.
    %   I   ... Scalar, denoting the source current.
    %   TX  ... Vector [2, x 1], denoting the source position.
    %
    % OUTPUT PARAMETER
    %   eleVS_fun ... Struct, containing function and gradient handles.

    %% Check input.

    assert(isscalar(rho), 'rho - Scalar denoting resistivity, expected.');
    assert(isscalar(I), 'I - Scalar denoting source current, expected.');
    assert(isvector(TX) && length(TX) == 2, ...
        'TX - Vector [2 x 1], denoting source position, expected.');

    %% Define function.

    eleVS_fun.f = @(X, Y) (rho * I) / (4 * pi * norm([X; Y] - TX(:)));

    %% Define gradient.

    if license('test', 'symbolic_toolbox')
        x_sym = sym('x', 'real');
        y_sym = sym('y', 'real');
        eleVS_fun.grad = [diff(eleVS_fun.f, x_sym); ...
                          diff(eleVS_fun.f, y_sym)];
        eleVS_fun.grad = matlabFunction(eleVS_fun.grad, ...
                            'Vars', {'x', 'y'});
        eleVS_fun.J = eleVS_fun.grad;
        clear('x_sym', 'y_sym');
    else
        % Skip derivation.
    end

    %% Set required quadrature order.

    eleVS_fun.quad_ord = 4;
end
