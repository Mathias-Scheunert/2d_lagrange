function DtN_fun = getD2NBCVals(TX)
    % Get kernel value handle for Dirichlet-to-Neumann operator.
    %
    % The idea is to formulate a mixed / Robin BC exploiting the known 
    % behavior of a electrode on top of a half-space.
    % 
    % u_ref(x,y)   = I / (pi * sigma) * K_0(k, r)
    % Grad_(u_ref) = - I / (pi * sigma) * k * K_1(k, r) * ((x, y) - TX) / r]
    % Grad_(u_ref) = - [((x, y) - TX) * k / r * K_1(k, r) / K_0(k, r)] * u_ref
    %
    % therefore
    %   d_(u_ref) / d_n                    = n * Grad_(u_ref)
    %   d_(u_ref) / d_n - n * Grad_(u_ref) = 0
    %
    % now one sets
    %   g_Robin = d_(u(x)) / d_(n) = [n * ((x,y) - TX)] / [log(r) * r^2] * u(x)
    % i.e. a Robin boundary condition without an additional Neumann part is
    % derived!
    %
    % with the kernel
    %   DtN_fun = ((x, y) - TX) * k / r * K_1(k, r) / K_0(k, r)
    %
    % with
    %   n     ... Vector of normal direction
    %   (x,y) ... Position vector
    %   TX    ... Position vector of source
    %   r     ... Distance (x,y) zu TX, d.h. ||(x,y) - TX||
    %
    % SYNTAX
    %   DtN_fun = getD2NVals(X, Y, TX)
    %
    % INPUT PARAMETER
    %   TX ... Vector [1 x 2] of source position.
    %   k  ... Scalar, wavenumber.
    %
    % OUTPUT PARAMETER
    %   DtN_fun ... Struct, containing function handle and quadrature order.          
    %
    % REMARKS
    % See Dey & Morrision, 1979
    %     Dissertation, Julia Weißflog, 2017, pp. 20
    %
    % TODO: The formulation only holds, if the TX is actually placed
    %       at the top of the surface!
    %       (for borehole measurement, see Bing, Z; 1998)

    % Check input.
    assert(isvector(TX), 'Only single source position supported.');
    
    % Set handles.
    r = @(x, y) sqrt((x - TX(1))^2 + (y - TX(2))^2);
    % Note: DtN_fun is set as an factory function, such that the parameter
    % k can be modified on demand.
    % I.e. the entire expression does not have to be fully defined 
    % in the DRIVE_ script but can rather be passed as a kind of template
    % to assembleDC25D.m and treatBC.m.
    % Therein, the wavenumber is finally set.
    DtN_fun = struct();
    DtN_fun.f = @(k) @(x, y) ([x; y] - TX(:)) * ...
                      k / r(x, y) * ...
                      getBesselKernel(k, r(x, y));

    % Set quadrature order.
    % TODO: fix - is that the correct order?
    DtN_fun.quad_ord = 4;
end

function kernel = getBesselKernel(k, r)
    % Evaluate ratio.
    kernel = besselk(1, k * r) / besselk(0, k * r);
    
    % Handle NaN.
    % TODO: fix - is that a hack or legit?
    if isnan(kernel)
        kernel = 1;
    end
end