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
    % TODO: check whether a sqrt(2) has to be added to [y_hat;z_hat].
    % -> Ensures the normal component of the function to have length 1,
    % like it is the case for the other two functions.
    % (but makes the transition to the other functions inconsistent at the
    %  corners of the edge).
    base.Phi = @(y_hat, z_hat) ...
                   [y_hat,     y_hat, y_hat - 1;
                    z_hat - 1, z_hat, z_hat];

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