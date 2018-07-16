function [err_L2, err_H1] = getError(mesh, fe, u, ref)
    % Calculate discretization error w.r.t. an analytic reference solution.
    %
    % The (square of the) error(-norms)s are defined by:
    % || u_exact - u_FE ||^2_L2 = L2 norm ^2 ...
    %   \int_omega (u_exact - u_FE)^T * (u_exact - u_FE) dx
    %
    % || grad(u_exact - u_FE) ||^2_L2 = H1-semi-norm ^2 = energy norm ...
    %   \int_omega \grad(u_exact - u_FE)^T * \grad(u_exact - u_FE) dx
    % 
    % || u_exact - u_FE ||_H1 = ...
    %   \sqrt(|| u_exact - u_FE ||^2_L2 + || \grad(u_exact - u_FE) ||^2_L2)
    %
    % SYNTAX
    %   [err_L2, err_H1] = getError(mesh, fe, u, ref)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   u    ... Vector, containing the solution of the FE-FWP, i.e. the
    %            FE-coefficients at the DOF
    %   ref  ... Struct, containing the information about the reference
    %            solution and its derivatives.
    %
    % OUTPUT PARAMETER
    %   err_L2 ... Scalar, denoting the discretization error 
    %              w.r.t L2-norm.
    %   err_H1 ... Scalar, denoting the discretization error 
    %              w.r.t H1-norm.
    %
    % REMARKS
    %   In order to provide the correct discretization error and not only
    %   calculating the error w.r.t. the interpolant (= point-wise
    %   comparison between the ref. sol and the FE sol. at the DOF) the
    %   integral norms have to be evaluated with an adapted quadrature
    %   order w.r.t the polynomial degree of the analytic reference
    %   solution.

    %% Check input
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isvector(u) && length(u) == fe.sizes.DOF, ...
        'u - FE solution vector expected.');
    assert(isstruct(ref) && all(isfield(ref, {'f', 'quad_ord'})), ...
        'ref - struct, including all information of reference solution, expected.');
    
    % Get common sizes.
    n_cell = fe.sizes.cell;
    
    %% Calculate discretization errors.
           
    % Get quadrature rule.
    [loc_quad_points, weights] = Quad.getQuadratureRule(ref.quad_ord, 2);
    sqrt_weights = num2cell(sqrt(weights));
        
    % Set up recurring quantity.
    % Get basis functions for all quadrature nodes referred to
    % the reference simplex.
    phi_quad = arrayfun(@(x,y) {fe.base.Phi(x, y)}, ...
        loc_quad_points(:,1), loc_quad_points(:,2));
    
    % Get basis function gradients for all quadrature nodes referred to
    % the reference simplex.
    grad_phi_quad = arrayfun(@(x,y) {fe.base.grad_Phi(x, y)}, ...
        loc_quad_points(:,1), loc_quad_points(:,2));
    
    % Iterate over all simplices.
    err_L2 = 0;
    err_H1_semi = 0;
    for ii = 1:n_cell     
       
        % Get global quadrature points.
        global_quad_points = (mesh.maps{ii}.B * loc_quad_points.' ...
                                + mesh.maps{ii}.b).';

        % Get reference solution at all quadrature nodes in global 
        % coordinates.
        u_ref = arrayfun(@(x, y) ref.f(x, y), ...
            global_quad_points(:, 1), global_quad_points(:, 2));
                
        % Get FE solution at all quadrature nodes in barycentric 
        % coordinates by superposing all basis function evaluations via
        % inner product.
        u_fe = cell2mat(cellfun(@(x) {x * u(fe.DOF_maps.cell2DOF{ii})}, ...
            phi_quad));
               
        % Calculate the square of the L2-norm of the error for the current 
        % simplex:
        % || u_fe - u_ref ||_L2^2 = \int_omega (u_fe - u_ref)^2 = ...
        %   \sum_k \abs(\det(B_k)) \sum_j weights_j ...
        %       (\sum_i u_i \phi_i(x_j) - u_ref(x_j))^2
        % k - num simplices
        % j - num quadrature nodes
        % i - num basis functions
        term1 = u_fe - u_ref;
        err_L2_loc = abs(mesh.maps{ii}.detB) * (weights .* term1).' * term1;
        
        % Add up all simplices.
        err_L2 = err_L2 + err_L2_loc;
        
        % Calculate square of H1-(seminorm) of the error for the current 
        % simplex:
        % || \grad (u_fe - u_ref) ||_L2^2 = \int_omega \grad(u_fe - u_ref)^2 = ...
        %   \sum_k \abs(\det(B_k)) \sum_j weights_j ...
        %       (\sum_i u_i \inv(B_k) \grad(\phi_i(x_j)) - \grad(u_ref(x_j)))^2
        % k - num simplices
        % j - num quadrature nodes
        % i - num basis functions
        if any(isfield(ref, {'J', 'grad'}))

            % Get gradient of reference solution at all quadrature nodes in 
            % global coordinates.
            if isfield(ref, 'grad')
                grad_u_ref = arrayfun(@(x, y) {ref.grad(x, y)}, ...
                    global_quad_points(:, 1), global_quad_points(:, 2));
            else
                grad_u_ref = arrayfun(@(x, y) {ref.J(x, y).'}, ...
                    global_quad_points(:, 1), global_quad_points(:, 2));
            end
            
            % Get FE solution gradient at all quadrature nodes in 
            % barycentric coordinates.
            grad_u_fe = cellfun(@(x) {...
                mesh.maps{ii}.BinvT * x * u(fe.DOF_maps.cell2DOF{ii})}, ...
                grad_phi_quad);

            % Calculate local error term.
            % As the gradient is nonscalar, the quadrature summation
            % is evaluated first (contrary to the L2-error handling).
            term1 = cellfun(@minus, grad_u_fe, grad_u_ref, ...
                'UniformOutput', false);
            term1 = cellfun(@(x, y) {x .* y}, term1, sqrt_weights);
            term1 = sum(cat(3, term1{:}), 3);
            err_H1_loc = abs(mesh.maps{ii}.detB) * (term1.' * term1);
            
            % Add up all simplices.
            err_H1_semi = err_H1_semi + err_H1_loc;
        else 
            err_H1_semi = 0;
        end
    end
    
    % Take square roots.
    err_H1 = sqrt(abs(err_H1_semi) + abs(err_L2));
    err_L2 = sqrt(abs(err_L2));
end

%{
function h = getDiam(mesh)
    % Calculate the max. inscribed circle diameter of simplices in mesh.
    %
    % r = 2*A / (a+b+c) = sqrt((s-a)(s-b)(s-c) / s) with s = (a+b+c) / 2
    % h = 2*r
    %
    % OUTPUT PARAMETER
    %   h      ... Scalar, denoting the biggest diameter of the inscribed 
    %              circles for all triangles/simplices in mesh.
    % TODO: implement or discard.
    
end
%}