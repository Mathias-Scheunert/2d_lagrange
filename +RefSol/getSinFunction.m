function sin_fun = getSinFunction()
    % Get information for 2D sin and polynomial function.
    %
    % Provides relevant information of an analytic sin-function:
    %   f(x,y) = −10x² + 10y² + 4sin(xy) − 2x + x⁴
    %
    % SYNTAX
    %   sin_fun = getSinFunction()
    %
    % OUTPUT PARAMETER
    %   sin_fun ... Struct, containing function handles of the sin function
    %               as well as its Jacobian and Hessian matrix.
    
    %% Define fuction.

    sin_fun.f = @(X, Y) -10 * X .^ 2 + 10 * Y .^ 2 + ...
        4 * sin(X .* Y) - 2 * X + X.^ 4;

    %% Define derivatives.

    % Define Jacobian terms.
    Jx = @(X, Y) -20 * X + 4 * Y .* cos(X .* Y) - 2 + 4 * X .^ 3;
    Jy = @(X, Y) 20 * Y + 4 * X .* cos(X .* Y);

    % Summarize Jacobian (row vector per definition).
    sin_fun.J = @(X, Y) [Jx(X, Y), Jy(X, Y)];

    % Define Hessian terms.
    Hxx = @(X, Y) -20 - 4 * Y .^ 2 .* sin(X .* Y) + 12 * X .^ 2; 
    Hxy = @(X, Y) 4 * cos(X .* Y) - 4 * X .* Y .* sin(X .* Y);
    Hyx = @(X, Y) 4 * cos(X .* Y) - 4 * X .* Y .* sin(X .* Y);
    Hyy = @(X, Y) 20 - 4 * X .^ 2 .* sin(X .* Y);

    % Summarize Hessian matrix and -Laplace.
    sin_fun.H = @(X, Y) [Hxx(X, Y), Hxy(X, Y); ...
                         Hyx(X, Y), Hyy(X, Y)];
    sin_fun.L = @(X, Y) -(Hxx(X, Y) + Hyy(X, Y));
    
    % Set required quadrature order.
    sin_fun.quad_ord = 6;
end
