function Peak = getPeakFunction(eps_peak, pos)
    % Get information for 2D peak function which simulates Dirac.
    %
    % SYNTAX
    %
    %   Peak = getPeakFunction(eps_peak)
    %
    % INPUT PARAMETER
    %   eps_peak ... Scalar, denoting the amplitude and the width of the
    %                normal distribution function which is used to 
    %                approximate the Dirac impulse.
    %   pos      ... Vector, providing the peak position.
    %
    % OUTPUT PARAMETER
    %
    %   Peak ... Struct, containing function handle.
    %
    % REMARKS
    %   See https://de.wikipedia.org/wiki/Mehrdimensionale_Normalverteilung
    
    %% Check input
    
    assert(isscalar(eps_peak), ...
        'eps_peak - scalar, denoting the amplitude, expected.');
    assert(isvector(pos) && length(pos) == 2, ...
        'pos - vector, denoting the peak position, expected.');
    
    %% Define peak-fuction.
    
    % Arbitrarily modify shape of the peak.
    % The peak should comprise a somehow small area and have the desired
    % amplitude.
    % TODO: automatically determine scaling factors w.r.t. the cell sizes.
    scale_C = eps_peak * 1e-4;
    scale_f = scale_C ^ 1.65;
    
    % Define function (assuming a uniform covariance matrix).
    cord = @(x, y) [x - pos(1); y - pos(2)];
    C = scale_C * eye(2, 2);
    CI = 1/scale_C * eye(2, 2);
    Peak.f = @(x, y) scale_f * 1/(2 * pi * det(C)) * ...
        exp(-(cord(x,y).' * CI * cord(x,y)) / 2);
end