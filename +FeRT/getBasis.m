function base = getBasis()
    % Provide Raviart-Thomas zeroth order vectorial basis.
    %
    % SYNTAX
    %   base = getBasis()
    %
    % OUTPUT PARAMETER
    %   base ... Struct, including the function handles to the basis 
    %            functions and their divergence.
        
    %% Get function and divergence.
    
    % Define basis functions.
    base.Phi = @(y_hat, z_hat) ...
                   [y_hat,    y_hat, y_hat - 1;
                    z_hat - 1,z_hat, z_hat];

    % Define basis function divergence.
    base.div_Phi = @(y_hat, z_hat) ...
                      [2, 2, 2];

    % Set local DOF coordinates for the basis functions.
    base.DOF = [0.5, 0;
                0.5, 0.5;
                0,   0.5];
    
    % Set local DOF direction (normal vector).
    base.DOF_orient = [ 0, -1;
                        1,  1;
                       -1,  0];
end