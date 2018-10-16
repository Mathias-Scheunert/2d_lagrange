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
        
    % Set up struct.
    base = struct();
    base.name = 'Raviart-Thomas';
    
    % Define basis functions.
    % TODO: check whether a sqrt(2) has to be added to [x_hat; y_hat].
    % -> I.g. Ensures the normal component of the function to have norm=1,
    % like it is the case for the other two functions.
    % (but makes the transition to the other functions inconsistent at the
    %  corners of the edge).
    % In practice however, all normal components are already similar?!
    base.Phi = @(x_hat, y_hat) ...
                   [x_hat,     x_hat, x_hat - 1;
                    y_hat - 1, y_hat, y_hat];
%     base.Phi = @(x_hat, y_hat) ...
%                    [x_hat,     sqrt(2) * x_hat, x_hat - 1;
%                     y_hat - 1, sqrt(2) * y_hat, y_hat];

    % Define basis function divergence.
    base.div_Phi = @(x_hat, y_hat) ...
                      [2, 2, 2];

    % Set local DOF coordinates for the basis functions.
    base.DOF = [0.5, 0;
                0.5, 0.5;
                0,   0.5];
    
    % Set local DOF direction (normal vector).
    base.DOF_normals = [[0, -1];
                        [1,  1] / sqrt(2);
                       [-1,  0]];
end