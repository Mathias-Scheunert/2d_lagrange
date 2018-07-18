function poisson_fun = getPoisson2D(TX)
    % Analytic solution of Poisson equation in 2D:
    %
    %   \u(r) = -ln(|r|)/(2 \pi)
    %         r = observation point = \sqrt(x^2 + y^2 + z^2)
    %
    % SYNTAX
    %   poisson_fun = getPoisson2D(TX)
    %
    % INPUT PARAMETER
    %   TX ... Vector [2, x 1], denoting the source position.
    %
    % OUTPUT PARAMETER
    %   poisson_fun ... Struct, containing function and gradient handles.
    
    %% Define function.

    poisson_fun.f = @(x, y) -1 / (2 * pi) * log(norm([x; y] - TX(:)));

    %% Define gradient.

    if license('test', 'symbolic_toolbox')
        x_sym = sym('x', 'real');
        y_sym = sym('y', 'real');
        poisson_fun.grad = [diff(poisson_fun.f, x_sym); ...
                            diff(poisson_fun.f, y_sym)];
        poisson_fun.grad = matlabFunction(poisson_fun.grad, ...
                              'Vars', {'x', 'y'});
        poisson_fun.J = poisson_fun.grad;
        clear('x_sym', 'y_sym');
    else
        % Skip derivation.
    end
    
    %% Set required quadrature order.
    
    poisson_fun.quad_ord = 4;
end