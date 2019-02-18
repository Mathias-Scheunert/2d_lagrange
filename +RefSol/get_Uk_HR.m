function Uk = get_Uk_HR(k, rTXRX, I, rho)
    % Calculates the line source 2D DC potential for a hom. half-space.
    % 
    % SYNTAX
    %   Uk = get_Uk_2layer_2D(k, rTXRX, h, I, rho1, rho2)
    %
    % INPUT PARAMETER
    %   rTXRX  ... Vector [k x 1] of TX-RX offsets (=|rTX - rRX|).
    %   h      ... Thickness of first layer.
    %   I      ... TX current / strengh.
    %   rho1   ... Specific resistivity of layer.
    %
    % OUTPUT PARAMETER
    %   Uk ... Matrix [length(rTXRX), length(k)] of solutions.

    % Fetch infos.
    n_r = length(rTXRX);
    n_k = length(k);
    
    % Calculate Uk.
    Uk = zeros(n_r, n_k);
    for jj = 1:n_k
        Uk(:,jj) = (I * rho) / (2*pi) * ...
                   arrayfun(@(r) besselk(0, k(jj) * r), rTXRX);
    end
end