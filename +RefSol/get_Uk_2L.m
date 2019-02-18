function Uk = get_Uk_2L(k, rTXRX, h, I, rho1, rho2)
    % Calculates the line source 2D DC potential for a 2-layered earth.
    % 
    % SYNTAX
    %   Uk = get_Uk_2L(k, rTXRX, h, I, rho1, rho2)
    %
    % INPUT PARAMETER
    %   rTXRX  ... Vector [k x 1] of TX-RX offsets (=|rTX - rRX|).
    %   h      ... Thickness of first layer.
    %   I      ... TX current / strengh.
    %   rho1/2 ... Specific resistivity of layers.
    %
    % OUTPUT PARAMETER
    %   Uk ... Matrix [length(rTXRX), length(k)] of solutions.
    %
    % COPYRIGHT
    %
    %   Function originally introduced by Z.Bing, 1998, Dissertation.
    
    % Fetch infos.
    n_r = length(rTXRX);
    n_k = length(k);
    
    % Set reflection coefficient.
    R12 = (rho2 - rho1) / (rho2 + rho1);
    
    % Calculate Uk.
    Uk = zeros(n_r, n_k);
    for ii = 1:n_k
        Uk(:,ii) = (I * rho1) / (2*pi) * (besselk(0, k(ii) * rTXRX) ...
                    + 2 * get_summation(h, k(ii), rTXRX, R12));
    end
end

function sum = get_summation(h, k, rTXRX, R12)
    % Approximate the infinite summation.

    % Set up problem.
    max_iter = 1000;
    sum = R12 * besselk(0, k * sqrt(rTXRX.^2 + (2*h)^2)); % initial summand
    tol = 1e-5 * norm(sum);
    
    % Loop until convergence is reached.
    for ii = 2:max_iter
       sum_add = R12^ii * besselk(0, k * sqrt(rTXRX.^2 + (2*ii*h)^2));
       if norm(sum_add) < tol
           break;
       end
       sum = sum + sum_add;
    end
end


