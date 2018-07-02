function peak_fun = getPeakFunction(val, pos, eps_peak)
    % Get information for 2D peak function which simulates Dirac.
    %
    % Arbitrarily modify shape of the peak.
    % The peak should comprise a somehow small area compared to h.
    %
    % SYNTAX
    %
    %   peak_fun = getPeakFunction(eps_peak)
    %
    % INPUT PARAMETER
    %   val      ... Scalar [or vector] scaling factor[s], denoting source 
    %                strength.
    %   pos      ... Vector [or matrix], providing the peak position[s].
    %   eps_peak ... Scalar, controlling the amplitude and the width of the
    %                normal distribution function which is used to 
    %                approximate the Dirac impulse.
    %
    % OUTPUT PARAMETER
    %
    %   peak_fun ... Struct, containing function handle.
    %
    % REMARKS
    %   See https://de.wikipedia.org/wiki/Mehrdimensionale_Normalverteilung
    %
    % TODO: automatically determine scaling factors w.r.t. the cell sizes.
    
    %% Check input
    
    assert(isscalar(val) || isvector(val), ...
        'val - scalar [or vector], denoting the amplitude[s], expected.');
    assert((isvector(pos) && length(pos) == 2) || ...
        (ismatrix(pos) && size(pos, 2) == 2), ...
        'pos - vector [or matrix], denoting the peak position[s], expected.');
    assert(isscalar(eps_peak), ...
        'eps_peak - scalar, controlling the width, expected.');
    assert(length(val) == size(pos, 1), ...
        'Length of val has to coincide with rows in pos.')
    
    %% Define peak-fuction.
        
    % Get sizes.
    n_peak = length(val);
    
    % Define function (assuming a uniform covariance matrix).
    C = eps_peak * eye(2, 2);
    CI = 1/eps_peak * eye(2, 2);
    peak_fun.f = @(x, y) 0;
    for ii = 1:n_peak
        cord = @(x, y) [x - pos(ii, 1); y - pos(ii, 2)];
        peak_fun.f = @(x, y) peak_fun.f(x,y) +  ...
            val(ii) * (1/(2 * pi * det(C)) * ...
            exp(-(cord(x,y).' * CI * cord(x,y)) / 2));
    end
    
    % Set required quadrature order.
    peak_fun.quad_ord = 3;
end