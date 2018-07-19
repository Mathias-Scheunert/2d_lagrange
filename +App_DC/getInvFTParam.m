function params = getInvFTParam(TX, RX, type)
    % Get weights and or nodes of the quadrature approach for inverse FT.
    %
    % SYNTAX
    %   params = getInvFTParam(TX, RX, type)
    %
    % INPUT PARAMETER
    %   TX   ... Vector [1 x 2], denoting the source position.
    %   RX   ... Matrix [n x 2], denoting the receiver positions.
    %   type ... Char, denoting the quadrature type.
    %            type = {'Boerner', 'Bing', 'Xu'}
    %
    % OUTPUT PARAMETER
    %   params ... Struct, containing the relevant quadrature approach
    %              information. See respective subfunctions for further
    %              details.
    
    %% Check input.
    
    assert(isvector(TX) && length(TX) == 2, ...
        'TX - Vector [1 x 2] denoting source position, expected.');
    assert(ismatrix(RX) && size(RX, 2) == 2, ...
        'RX - Matrix [n x 2] denoting receiver positions, expected.');
    assert(ischar(type), ...
        'type - Char, denoting the quadrature type, expected.')
    
    %% Get infos.
    params = struct();
    switch type
        case 'Boerner'
            [params.k, params.w, params.n] = getInvFTBoerner(TX, RX);
        case 'Bing'
            [params.k, params.n] = getInvFTBing(TX, RX);
        case 'Xu'
            [params.k, params.w, params.n] = getInvFTXu();
        otherwise
            error('Unknown type - "Boerner", "Bing" or "Xu" supported.')
    end
    params.type = type;
end

function [k, w, n] = getInvFTBoerner(TX, RX)
    % Calculate nodes and weights for quadrature approach by Boerner R.-U.
    %
    % Approach exploits logarithmic and exponetial asymptotic behaviour of
    % Bessel functions K_0(k * r), denoting the analytic solution for a 
    % unit point source over unit half-space.
    %
    % Note that source strengh I and resistivity \rho only shift the
    % amplitude of the analytic solution such that dealing with unit
    % quantities is justified.
    %
    % SYNTAX
    %   [k, w, n] = getInvFTBoerner(TX, RX)
    %
    % OUTPUT PARAMETER
    %   k ... Vector, containing quadrature nodes (spatial wavenumbers).
    %   w ... Vector, containing weights for quadrature summation.
    %   n ... Scaler, denoting number of k or w, respectively.
    
    
    %% Prepare parameter.
    
    % Get TX-RX offsets.
    r_TX2RX = sqrt((TX(1) - RX(:,1)).^2 + (TX(2) - RX(:,2)).^2);
    
    % Get critical wavenumber (transition of asymptotic behavior).
    k_crit = 1 / (2 * min(r_TX2RX));

    % Define number of wavenumbers for both parts of the asymptotic 
    % behavior of the analytic solution.
    % TODO: find usefull approach to somehow predetermine numbers of k's.
    n_k_log = 10;
    n_k_exp = 10;
    n = n_k_log + n_k_exp;

    % Get quadrature rules for both asymptotics.
    [x_Leg, w_Leg] = Quad.getQuadratureRule(n_k_log, 1, ...
        'type', 'Legendre', 'bnd', [0, 1]);
    [x_Lag, w_Lag] = Quad.getQuadratureRule(n_k_exp, 1, ...
        'type', 'Laguerre');

    % Derive wavenumbers and weights.
    k_log = k_crit * x_Leg(:).^2;
    w_log = 2 * k_crit * (x_Leg(:) .* w_Leg(:));
    k_exp = k_crit * (x_Lag(:) + 1);
    w_exp = k_crit * (exp(x_Lag(:)) .* w_Lag(:));
    k = [k_log; k_exp];
    w = [w_log; w_exp];
end

function [k, n] = getInvFTBing(TX, RX)
    % Calculate nodes and weights for quadrature approach by Bing Z. 1998.
    %
    % The values are based on the consideration of the bessel function
    % behaviour for variing wavenumbers w.r.t. the minimal and maximal
    % TX-RX offset.
    %
    % SYNTAX
    %   [k, n] = getInvFTBing(TX, RX)
    %
    % OUTPUT PARAMETER
    %   k   ... Vector, containing wavenumbers.
    %   n   ... Scalar, denoting number of k or w, respectively.
    
    debugging = pick(1, false, true);
    
    %% Prepare parameter.
    
    % Get TX-RX offsets.
    r_TX2RX = sqrt((TX(1) - RX(:,1)).^2 + (TX(2) - RX(:,2)).^2);
    r_min = min(r_TX2RX);
    r_max = max(r_TX2RX);

    % Get critical wavenumber (transition of asymptotic behavior).
    k_crit = 1 / (2 * r_min);
    
    % Set tolerance for curvature criteria.
    tol_min = 1e-3;
    tol_max = 1e-30;
    
    % Define number od wavenubers.
    n = 30;
        
    %% Get k range.
    
    % Predefine a large range of wavenumbers.
    k_range = logspace(-5, 4, 100);

    % Define analytic solution for unit source over unit half-space.
    I = 1;
    rho = 2 * pi;
    u_HS = @(k, r) (I * rho) / (2 * pi) * besselk(0, k * r);

    % Get solutions spectra and curvatures.
%     u_min = arrayfun(@(k) u_HS(k, r_min), k_range);
%     curv_u_min = gradient(gradient(u_min));
    u_max = arrayfun(@(k) u_HS(k, r_max), k_range);
    curv_u_max = gradient(gradient(u_max));
    
    % Get respective wavenumber boundaries.
    k_min = k_range(find(curv_u_max > tol_min, 1, 'first'));
%     k_max = k_range(find(curv_u_min > tol_max, 1, 'last'));
    
    % Get range of k and index for critical k.
%     k = logspace(log10(k_min), log10(k_max), n);
%     idx = find(k < k_crit, 1, 'last');
    k = logspace(log10(k_min), log10(k_crit), n);
    
    if debugging
        % Plot solutions.
        figure();
        loglog(k_range, u_min, 'b', k_range, u_max, 'r');
        hold on
        plot(k_min, ...
            u_max(find(curv_u_max > tol_min, 1, 'first')), 'x')
        plot(k_max, ...
            u_min(find(curv_u_min > tol_max, 1, 'last')), 'o')
        hold off
        legend(sprintf('r = %d', r_min), ...
               sprintf('r = %d', r_max), ...
               'k_{min}', 'k_{max}');
        ylim([1e-100, 1e5]);
        xlabel('k');
        ylabel('K_0(k, r)');
    end
end

function [k, w, n] = getInvFTXu()
    % Calculate nodes and weights for quadrature approach by Xu S.-Z. 2000.
    %
    % Nodes and weights are directly imported from the paper.
    %
    % SYNTAX
    %   [k, w, n] = getInvFTXu()
    %
    % OUTPUT PARAMETER
    %   k ... Vector, containing quadrature nodes (spatial wavenumbers).
    %   w ... Vector, containing weights for quadrature summation.
    %   n ... Scaler, denoting number of k or w, respectively.
    
    k = [0.0217102, 0.2161121, 1.0608400, 5.0765870].';
    w = [0.0463660, 0.2365931, 1.0382080, 5.3648010].';
    n = length(k);
end