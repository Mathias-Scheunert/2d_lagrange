function eleHS_fun = getElectrodeAtVS(rho, I, TX)
    % Get information for an point source in a homogeneous full-space.
    %
    % Provides relevant information of the potential of a point source 
    % in a homogeneous full-space in 3D:
    %   \Phi(r) = (\rho I) / (4 \pi \abs(r - r'))
    %         r = observation point = \sqrt(x^2 + y^2 + z^2)
    %         r'= source position
    %         I = source current
    %      \rho = resistivity of the half-space
    %
    % SYNTAX
    %   eleHS_fun = getElectrodeAtHS()
    %
    % OUTPUT PARAMETER
    %   eleHS_fun ... Struct, containing function handles of the sin function
    %                 as well as its Jacobian and Hessian matrix.
    
    %% Define fuction.
    
    pi = 3.141592653589793;
    r = @(X, Y) norm([X; Y] - TX(:));
    eleHS_fun.f = @(X, Y) (rho * I) / (4 * pi * r(X, Y));

    %% Define gradient.
    
    if license('test', 'symbolic_toolbox')
        x_sym = sym('x', 'real');
        y_sym = sym('y', 'real');
        eleHS_fun.grad = [diff(eleHS_fun.f, x_sym); diff(eleHS_fun.f, y_sym)];
        eleHS_fun.grad = matlabFunction(eleHS_fun.grad, 'Vars', {'x', 'y'});
        eleHS_fun.J = eleHS_fun.grad;
        clear('x_sym', 'y_sym');
    else
        % Skip derivation.
    end

    % Set required quadrature order.
    eleHS_fun.quad_ord = 4;
end
