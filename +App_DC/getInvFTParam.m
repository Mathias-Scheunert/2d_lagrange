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
        ['TX - Vector [1 x 2] denoting source position, expected. ', ...
        'This is required in order to derive propper spatial wave ', ...
        'number selection']);
    assert(ismatrix(RX) && size(RX, 2) == 2, ...
        ['RX - Matrix [n x 2] denoting receiver positions, expected. ', ...
        'This is required in order to derive propper spatial wave ', ...
        'number selection']);
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
            [params.k, params.w, params.n] = getInvFTXu(TX, RX);
%             [params.k, params.w, params.n] = getInvFTXu_fix();
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
    n_k_log = 30;
    n_k_exp = 30;
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
    
    %% Prepare parameter.
    
    % Get TX-RX offsets.
    r_TX2RX = sqrt((TX(1) - RX(:,1)).^2 + (TX(2) - RX(:,2)).^2);
    r_min = min(r_TX2RX);
    r_max = max(r_TX2RX);

    % Get critical wavenumber (transition of asymptotic behavior).
    k_crit = 1 / (2 * r_min);
    
    % Set tolerance for curvature criteria.
    tol_min = 2e-1;
    
    % Define number od wavenubers.
    n = pick(2, 20, 60);
        
    %% Get k range.
    
    % Predefine a large range of wavenumbers.
    k_range = logspace(-5, 5, 100);

    % Define analytic solution for unit source over unit half-space.
    I = 1;
    rho = 2 * pi;
    u_HS = @(k, r) (I * rho) / (2 * pi) * besselk(0, k * r);

    % Get solutions spectra and curvatures.
    u_max = arrayfun(@(k) u_HS(k, r_max), k_range);
    grad_u_max = gradient(u_max);
    
    % Get respective wavenumber boundaries.
    k_min = k_range(find(abs(grad_u_max) < tol_min, 1, 'first'));
    k = logspace(log10(k_min), log10(k_crit), n);
end

function [k, w, n] = getInvFTXu(TX, RX)
    % Calculate weights by Xu S. 2000 for quadrature nodes by Bing Z 1998.
    %
    % The values are based on the consideration of the bessel function
    % behaviour for variing wavenumbers w.r.t. the minimal and maximal
    % TX-RX offset.
    %
    % SYNTAX
    %   [k, n] = getInvFTXu(TX, RX)
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
    
    % Set tolerance for curvature criteria.
    tol_min = 2e-1;
    tol_max = 1e-20;
    
    % Define number od wavenubers.
    n = pick(2, 20, 60);
        
    %% Get k range.
    
    % Predefine a large range of wavenumbers.
    k_range = logspace(-5, 5, 100);

    % Define analytic solution for unit source over unit half-space.
    I = 1;
    rho = 2 * pi;
    u_HS = @(k, r) (I * rho) / (2 * pi) * besselk(0, k * r);

    % Get solutions spectra and curvatures.
    u_min = arrayfun(@(k) u_HS(k, r_min), k_range);
    grad_u_min = gradient(u_min);
    u_max = arrayfun(@(k) u_HS(k, r_max), k_range);
    grad_u_max = gradient(u_max);
    
    % Get respective wavenumber boundaries.
    k_min = k_range(find(abs(grad_u_max) < tol_min, 1, 'first'));
    k_max = k_range(find(abs(grad_u_min) > tol_max, 1, 'last'));
    
    % Get range of k and index for critical k.
    k = (logspace(log10(k_min), log10(k_max), n)).';

    %% Get w and k by nonlinear least-squares.
%{ 
    % For a point source over a homogeneous half-space the solution
    %   (\rho * I) / (2 * pi * abs(r))
    % is ~ to
    %   (2 / pi) \int_{0}^{\inf} (I * rho) / (2 * pi) * besselk(0, k * r)
    % Hence, by general quadrature approach 
    %   1 / abs(r) ~ \sum_{j = 1}^{n} besselk(0, k_j * r) w_j
    % By defining 
    %   A*w = v
    % with
    %   A = [r_i besselk(0, k_j * r_i)]_{i=1:m,j:n}
    %   w = [w_j]_{j=1:n}
    %   v = [1]_{i=1:m}
    % a nonlinear minimization can be obtained to estimate optimal w_j and 
    % k_j by:
    %   ||v - A(k)*w||_2^2 -> min_{w,k}
    %   ||v - A(k_0)*w_0 - S_0 * \delta [k;w]||_2^2 -> min_{[k;w]}
    % with
    %   S(k_0,w_0)_0 = d_{A(k_0)*w_0} / d_{[k;w]}
    %          [k;w] = [k_0;w_0] + \delta [k;w]
    %   
    % Procedure: 
    %   1) guess k_0 -> Bing-approach and set w_0
    %   2) iterate until ||v - A(k)*w|| converges
    
    % Set up optimization for m = [k; w].
    lambda = 1e-5;
    m = log([k; 1/n * ones(n, 1)]);
    iter = 1;
    iter_max = 20;
    res = zeros(iter_max, 1);

    while iter <= iter_max
        
        % Get derivative w.r.t. log(w): (equals A since dependencies are 
        % linear)
        m_op = exp(m); % used within the subroutines.
        S_w = getSys(u_HS, m_op(1:n), r_TX2RX);
        
        % Get current residuum.
        res(iter) = norm(ones(length(r_TX2RX), 1) - S_w * m(n + 1:end));
        
        % Terminate if already converged.
        if (iter > 1) && (abs(res(iter) - res(iter - 1)) < abs(res(1))*1e-3)
            break;
        end
        
        % Get derivative w.r.t. log(k): (using perturbation method)
        S_k = getSens(u_HS, m_op(1:n), r_TX2RX, ...
                  S_w * m(n + 1:end), m(n + 1:end));
        % (include diagonal matrix, resulting from chain rule)
        S_k = S_k * diag(exp(m(1:n)));
    
        % Set up linearized system for optimization of k and w.
        S = [S_k, S_w];
        b = (1 + zeros(length(r_TX2RX), 1)) - (S_w * m(n + 1:end));

        % Solve corresponding normal equation (Levenberg-Maquard-scheme).
        dm =  ((S.'*S + lambda * eye(2*n, 2*n)) \ (S.' * b));
        m = m + dm;
        iter = iter + 1;
    end
    
    % Exclude wavenumbers and weights.
    % (Add prefactor to derive weights, equivalent to Boerner and Bing.)
    m = exp(m);
    k = m(1:n);
    w = (pi / 2) * m(n + 1:end);
%}   
    %% Get corresponding w by linear least-squares.
% %{
    % For a point source over a homogeneous half-space the solution
    %   (\rho * I) / (2 * pi * abs(r))
    % is ~ to
    %   (2 / pi) \int_{0}^{\inf} (I * rho) / (2 * pi) * besselk(0, k * r)
    % Hence, by general quadrature approach 
    %   1 / abs(r) ~ \sum_{j = 1}^{n} besselk(0, k_j * r) w_j
    % By defining 
    %   A*w = v
    % with
    %   A = [r_i besselk(0, k_j * r_i)]_{i=1:m,j:n}
    %   w = [w_j]_{j=1:n}
    %   v = [1]_{i=1:m}
    % a linear minimization can be obtained to estimate optimal w_j for a
    % set of given k_j by:
    %   ||v - A*w||_2^2 -> min_w
    %   w = \inv(A^{T}A) * A^{T}*v
    %
    % Procedure: 
    %   1) take k -> Bing-approach
    %   2) linear optimize for w
    
    lambda = 1e-6;
    
    % Set up linear system for estimation of weights.
    A_w = getSys(u_HS, k, r_TX2RX);
    b_w = 1 + zeros(length(r_TX2RX), 1);

    % Solve corresponding normal equation (Levenberg-Maquard-scheme).
    % (Add prefactor to derive weights,  equivalent to Boerner and Bing.)
    w = (pi / 2) * ((A_w.'*A_w + lambda * eye(n, n)) \ (A_w.' * b_w));
% %}

    %% Plot asymptotics.
    
    if debugging
        figure();
        loglog(k_range, u_min, 'b', k_range, u_max, 'r');
        hold on
        plot(k_min, ...
            u_max(find(abs(grad_u_max) < tol_min, 1, 'first')), 'x')
        plot(k_max, ...
            u_min(find(abs(grad_u_min) > tol_max, 1, 'last')), 'o')
        hold off
        legend(sprintf('r = %d', r_min), ...
               sprintf('r = %d', r_max), ...
               'k_{min}', 'k_{max}', ...
               'Location', 'best');
        ylim([1e-100, 1e5]);
        xlabel('k');
        ylabel('K_0(k, r)');
    end
end

function A = getSys(fun, k, r)
    % Calculates system matrix of linear problem w.r.t. w.
   
    A = zeros(length(r), length(k));
    for ii = 1:length(r)
        A(ii,:) = abs(r(ii)) * arrayfun(@(x) fun(x,r(ii)), k);
    end
end

function S = getSens(fun, k, r, Aw_org, w)
    % Calculates d_(A*w) / d_k by perturbation.
    %
    % TODO: add Taylor test for verification.
    
    % Initialize.
    S = ones(length(r), length(k));
    noise = 1e-6;

    % Loop over columns.
    for ii = 1:length(k)
        
        % Set current k_dist to original and add noise.
        k_dist = k;
        k_dist(ii) = k(ii).*(1 + noise);
        
        % Get corresponding solution.
        A_dist = getSys(fun, k_dist, r);
        Aw_dist = A_dist * w;
        
        % Calculate colums of sensitivity.
        S(:,ii)= (Aw_dist - Aw_org) / (k_dist(ii) * noise);
    end
end