function eleHS_fun = getElectrodeAtHS(rho, I, TX)
    % Get information for an point source at a homogeneous half-space.
    %
    % Provides relevant information of the potential of a point source 
    % above a homogeneous half-sapce in 3D:
    %   \Phi(r) = (\rho I) / (2 \pi \abs(r - r'))
    %         r = observation point = \sqrt(x^2 + y^2 + z^2)
    %         r'= source position
    %
    % SYNTAX
    %   eleHS_fun = getElectrodeAtHS(rho, I, TX)
    %
    % INPUT PARAMETER
    %   rho ... Scalar, denoting the HS resistivity. 
    %   I   ... Scalar, denoting the source current.
    %   TX  ... Vector [2, x 1], denoting the source position.
    %
    % OUTPUT PARAMETER
    %   eleHS_fun ... Struct, containing function and gradient handles.
    
    %% Check input.
    
    assert(isscalar(rho), 'rho - Scalar denoting resistivity, expected.');
    assert(isscalar(I), 'I - Scalar denoting source current, expected.');
    assert(isvector(TX) && length(TX) == 2, ...
        'TX - Vector [2 x 1], denoting source position, expected.');
    
    %% Define function.
    
    r = @(X, Y) norm([X; Y] - TX(:));
    eleHS_fun.f = @(X, Y) (rho * I) / (2 * pi * r(X, Y));

    %% Define gradient.
    
    if license('test', 'symbolic_toolbox')
        x_sym = sym('x', 'real');
        y_sym = sym('y', 'real');
        eleHS_fun.grad = [diff(eleHS_fun.f, x_sym); ...
                          diff(eleHS_fun.f, y_sym)];
        eleHS_fun.grad = matlabFunction(eleHS_fun.grad, ...
                            'Vars', {'x', 'y'});
        eleHS_fun.J = eleHS_fun.grad;
        clear('x_sym', 'y_sym');
    else
        % Skip derivation.
    end

    %% Set required quadrature order.
    
    eleHS_fun.quad_ord = 4;
end
