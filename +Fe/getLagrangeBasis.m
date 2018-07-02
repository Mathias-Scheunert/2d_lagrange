function base = getLagrangeBasis(order) 
    % Provide Lagrange basis coefficients using Vandermonde matrix.
    %
    % Point coordinates are referred to the two barycentric coordinates.
    % I.e. the kartesian coordinates of the 2D reference triangle /
    % simplex with 
    %   x_hat = [0, 1], y_hat = [0, 1] 
    % and the respective vertices 
    %   p1 = (0,0), p2 = (1,0), p3 = (0,1).
    %
    % For the Vandermonde matrix V and a polynomial basis functions it 
    % holds that
    %   V * [a_1, ... a_n].' = kron_i
    %
    % 1st order basis function, linear (in 2D):
    %   a_0 + a_1 * x + a_2 * y
    % -> 3 unknowns = 3 DOF at the element (triangle)
    %
    % 2nd order basis function, quadratic (in 2D):
    %   a_0 + a_1 * x + a_2 * y + a_3 * x^2 + a_4 * y^2
    % -> 6 unknowns = 6 DOF at the element
    %
    % SYNTAX
    %   base_coef = LagrangeBasis(order) 
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
    
    switch order
        case 1
            % First order Vandermonde matrix.
            % V * [a_1, ... a_n].' = kron_i
            V1 = [1 0 0;
                  1 1 0;
                  1 0 1];

            % Solving V * [a_1, ... a_n].' = kron_i for all i, provides 
            % the i basis function coefficients:
            base.coef = V1 \ eye(size(V1));

            % Define basis functions.
            Phi_var = @(x_hat, y_hat) ...
                [1; x_hat; y_hat];
            base.Phi = @(x_hat, y_hat) ...
                sum(bsxfun(@times, base.coef, Phi_var(x_hat, y_hat)));
    
            % Define basis function gradient.
            base.grad_Phi = @(x_hat, y_hat) reshape([...
                -1.0, -1.0, 1.0, 0.0, 0.0, 1.0], ...
                [2, 3]);
        
            % Get local DOF coordinates for the basis functions.
            base.DOF = V1(:,2:3);
            
            % Therefore derive the association between the index of the
            % local/reference simplex vertices to the global/mapped ones.
            % (See Mesh.getAffineMap for the definitions of Bk and bk.)
            base.DOF_loc2DOF_glo = [3, 1, 2];
             
        case 2
            % Second order Vandermonde matrix.
            V2 = [1 0   0   0    0    0;
                  1 1   0   1    0    0;
                  1 0   1   0    0    1;
                  1 0.5 0   0.25 0    0;
                  1 0.5 0.5 0.25 0.25 0.25;
                  1 0   0.5 0    0    0.25];
            base.coef = V2 \ eye(size(V2));

            % Define basis functions.
            Phi_var = @(x_hat, y_hat) ...
                [1; x_hat; y_hat; x_hat^2; x_hat * y_hat; y_hat^2];
            base.Phi = @(x_hat, y_hat) ...
                sum(bsxfun(@times, base.coef, Phi_var(x_hat, y_hat)));

            % Define basis function gradient. 
            base.grad_Phi =  @(x_hat, y_hat) reshape([...
                x_hat .* 4.0 + y_hat .* 4.0 - 3.0, ...
                x_hat .* 4.0 + y_hat .* 4.0 - 3.0, ...
                x_hat .* 4.0 - 1.0, ...
                0.0, ...
                0.0, ...
                y_hat .* 4.0 - 1.0, ...
                x_hat .* -8.0 - y_hat .* 4.0 + 4.0, ...
                x_hat .* -4.0, ...
                y_hat .* 4.0, ...
                x_hat .* 4.0, y_hat .* -4.0, ...
                x_hat .* -4.0 - y_hat .* 8.0 + 4.0], ...
                [2,6]);
            
            % Get local DOF coordinates for the basis functions.
            base.DOF = V2(:,2:3);
            
            % Therefore derive the association between the index of the
            % local/reference simplex vertices to the global/mapped ones.
            % See: 
            % Mesh.getAffineMap for the definitions of Bk and bk, referring
            %   to the nodes-DOF
            % Mesh.appendElementInfo for derivation of the edge index from
            %   the nodes
            % As the edge midpoints are not represented in the original
            % mesh structure, the edge index is used for that 
            % (See: Fe.getDOFMap).
            % Hence, the local DOF does not count from 1:6 but is rather
            % splitted in two counts from 1:3.
            % As the both, ordering of basis functions on local/reference 
            % simplex and ordering of edges on arbitrary triangle in mesh
            % follows the same principle, both vectors coincide.
            %              ind_node, ind_edge
            base.DOF_loc2DOF_glo = [[3, 1, 2], [3, 1, 2]];
    end
end