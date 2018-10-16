function base = getBasis(order) 
    % Provide Lagrange fist or second oder scalar basis.
    %
    % To derive expressions, Vandermonde matrix is exploited.
    %
    % Global point coordinates (x,y) are referred to the two barycentric 
    % coordinates (y,z)_hat (local coordinates). 
    % The third barycentric coordinate is linked to those two by
    % x_hat = 1 - y_hat - z_hat.
    % I.e. they can be treated as the kartesian coordinates of the 
    % 2D reference triangle / simplex with coordinate directions
    %   y_hat = [(0), 1, 0], z_hat = [(0), 0, 1] 
    % and the respective vertices (in "local kartesian coordinates") 
    %   p_local = p_local(y_hat, z_hat)
    %   p1 = (0,0), p2 = (1,0), p3 = (0,1).
    %
    % For the Vandermonde matrix V and a polynomial basis functions it 
    % holds that
    %   V * [a_1, ... a_n].' = kron_i
    %
    % 1st order basis function, linear (in 2D):
    %   a_0 + a_1 * y + a_2 * z
    % -> 3 unknowns = 3 DOF at the element (triangle)
    %
    % 2nd order basis function, quadratic (in 2D):
    %   a_0 + a_1 * y + a_2 * z + a_3 * y^2 + a_4 * z^2
    % -> 6 unknowns = 6 DOF at the element
    %
    % SYNTAX
    %   base = getBasis(order) 
    %
    % INPUT PARAMETER
    %   order ... Scalar, denoting the order of elements.
    %
    % OUTPUT PARAMETER
    %   base ... Struct, including all coefficients for the desired
    %            basis functions, the function handles to those and their
    %            gradients.
    %
    % REMARKS
    %   Definitions from 
    %       http://femwiki.wikidot.com/elements:lagrange-elements
     
    %% Check input.
    
    assert(isscalar(order) && order <= 2, ...
        'order - scalar 1 or 2 for first or second order linear basis functions.');
    
    %% Get Coefficients, functions, gradients.
    
    % Set up struct.
    base = struct();
    base.name = 'Lagrange';
    
    switch order
        case 1
            % First order Vandermonde matrix.
            % V * [a_1, ... a_n].' = kron_i
            V1 = [1 0 0;
                  1 1 0;
                  1 0 1];

            % Solving V * [a_1, ... a_n].' = kron_i for all i, provides 
            % the i basis function coefficients:
            coef = V1 \ eye(size(V1));

            % Define basis functions.
            Phi_var = @(y_hat, z_hat) ...
                [1; y_hat; z_hat];
            base.Phi = @(y_hat, z_hat) ...
                sum(bsxfun(@times, coef, Phi_var(y_hat, z_hat)));
    
            % Define basis function gradient.
            base.grad_Phi = @(y_hat, z_hat) reshape([...
                -1.0, -1.0, 1.0, 0.0, 0.0, 1.0], ...
                [2, 3]);
        
            % Get local DOF coordinates for the basis functions.
            base.DOF = V1(:,2:3);
                         
        case 2
            % Second order Vandermonde matrix.
            V2 = [1 0   0   0    0    0;
                  1 1   0   1    0    0;
                  1 0   1   0    0    1;
                  1 0.5 0   0.25 0    0;
                  1 0.5 0.5 0.25 0.25 0.25;
                  1 0   0.5 0    0    0.25];
            coef = V2 \ eye(size(V2));

            % Define basis functions.
            Phi_var = @(y_hat, z_hat) ...
                [1; y_hat; z_hat; y_hat^2; y_hat * z_hat; z_hat^2];
            base.Phi = @(y_hat, z_hat) ...
                sum(bsxfun(@times, coef, Phi_var(y_hat, z_hat)));

            % Define basis function gradient. 
            base.grad_Phi =  @(y_hat, z_hat) reshape([...
                y_hat .* 4.0 + z_hat .* 4.0 - 3.0, ...
                y_hat .* 4.0 + z_hat .* 4.0 - 3.0, ...
                y_hat .* 4.0 - 1.0, ...
                0.0, ...
                0.0, ...
                z_hat .* 4.0 - 1.0, ...
                y_hat .* -8.0 - z_hat .* 4.0 + 4.0, ...
                y_hat .* -4.0, ...
                z_hat .* 4.0, ...
                y_hat .* 4.0, z_hat .* -4.0, ...
                y_hat .* -4.0 - z_hat .* 8.0 + 4.0], ...
                [2,6]);
            
            % Get local DOF coordinates for the basis functions.
            base.DOF = V2(:,2:3);
    end
end