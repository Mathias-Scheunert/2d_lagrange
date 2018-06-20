function [err_L2, err_H1, h] = getError(mesh, fe, u, ref)
    % Calculate discretization error w.r.t. an analytic reference solution.
    %
    % The error(-norms)s are defined by:
    % || u_exact - u_FE ||_L2^2 = ...
    %   int_omega (u_exact - u_FE)^T*(u_exact - u_FE) dx
    %
    % || u_exact - u_FE ||_H1semi^2 = ...
    %   int_omega grad(u_exact - u_FE)^T*grad(u_exact - u_FE)dx
    % 
    % || u_exact - u_FE ||_H1 = ...
    %   sqrt(|| u_exact - u_FE ||_L2^2 + || u_exact - u_FE ||_H1semi^2)
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
    %   err_L2 ... Scalar, denoting the discretization error w.r.t
    %              L2-norm.
    %   err_H1 ... Scalar, denoting the discretization error w.r.t
    %              H1-norm.
    %   h      ... Scalar, denoting the biggest diameter of the inscribed 
    %              circles for all triangles/simplices in mesh.
    %
    % REMARKS
    %   In order to provide the correct discretization error and not only
    %   calculating the error w.r.t. the interpolant (= point-wise
    %   comparison between the ref. sol and the FE sol. at the DOF) the
    %   integral norms have to be evaluated with an adapted quadrature
    %   order w.r.t the polynomial degree of the analytic reference
    %   solution.
    %
    % TODO: implement.

    %% Check input
    
    %% Calculate L2-norm of the discretization error.
    
    %% Calculate H1-seminorm of the discretization error.
    
    %% Calculate H1-norm of the discretization error.
    
end

function h = getDiam(mesh)
    % Calculate the max. inscribed circle diameter of simplices in mesh.
    %
    % r = 2*A / (a+b+c) = sqrt((s-a)(s-b)(s-c) / s) with s = (a+b+c) / 2
    % h = 2*r
    
    % TODO: implement or discard.
    
end