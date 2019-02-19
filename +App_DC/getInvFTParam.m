function params = getInvFTParam(TX, RX, type)
    % Get weights and or nodes of the quadrature approach for inverse FT.
    %
    % For multiple TX a general set of weights w and wavenumbers k is
    % determined.
    % (Note: w and k can be determined for each m TX positions but this 
    % will lead to m unique linear systems.)
    %
    % SYNTAX
    %   params = getInvFTParam(TX, RX, type)
    %
    % INPUT PARAMETER
    %   TX   ... Vector [m x 2], denoting the source positions.
    %   RX   ... Matrix [n x 2], denoting the receiver positions.
    %   type ... Char, denoting the quadrature type.
    %            type = {'Boerner', 'Bing', 'Xu'}
    %
    % OUTPUT PARAMETER
    %   params ... Struct, containing the relevant quadrature approach
    %              information. See respective subfunctions for further
    %              details.
    %
    % REMARKS
    %   The 2.5D approach uses the bahavior of a known analytic solution
    %   (e.g. point source in 3D half-space) to set up the numerical 
    %   integration of the seperate 2D solutions (including choice of 
    %   wavenumbers and weights) which forms the 2.5D solution.
    %   This can only act as an approximation for the solution of an 
    %   arbitrary shaped underground!
    %   -> I.e. the solutions should be veryfied by comparing them to a 
    %   full 3D simulation.
    %   Furthermore, this implies that only one single point source can be 
    %   treated by this approach such that multi-pole arrangements have to
    %   be simulated by adding up solutions for several sources.
    %
    %   2.5D approach derived from:    Dey A., Morrison H.F.; 1979
    %   num. integration derived from: Ralph-Uwe Börner (per.communication)
    %                                  Bing Z.; 1998 (Dissertation)
    %                                  Xu S.; 2000
    %
    % TODO: Debug Bing-variant (where to start/stop summation?).
    
    %% Check input.
    
    assert(ismatrix(TX) && size(TX, 2) == 2, ...
        ['TX - Matrix [m x 2] denoting source position, expected. ', ...
        'This is required in order to derive propper spatial wave ', ...
        'number selection']);
    assert(ismatrix(RX) && size(RX, 2) == 2, ...
        ['RX - Matrix [n x 2] denoting receiver positions, expected. ', ...
        'This is required in order to derive propper spatial wave ', ...
        'number selection']);
    assert(ischar(type), ...
        'type - Char, denoting the quadrature type, expected.')
    
    %% Get infos.
    
    % Define number of wavenubers.
    % (Bing / Xu)
    n_k = 17; 
    
    % Iterate over TX positions.
    params = struct();
    switch type
        case 'Boerner'
            [params.k, params.w, params.n] = getInvFTBoerner(TX, RX);
        case 'Bing'
            error('Approach not fully proved - do not use!');
            [params.k, params.n] = getInvFTBing(TX, RX, n_k); %#ok
        case 'Xu'
            error('Approach misleading - do not use!');
            [params.k, params.w, params.n] = getInvFTXu(TX, RX, n_k); %#ok
        otherwise
            error('Unknown type - "Boerner", "Bing" or "Xu" supported.')
    end
    params.type = type;
end

function [k, w, n] = getInvFTBoerner(TX, RX)
    % Calculate nodes and weights for quadrature approach by Boerner R.-U.
    %
    % Approach exploits logarithmic and exponetial asymptotic behaviour of
    % Bessel functions K_0(u) with u = k * r denoting the analytic solution
    % for a unit point source over unit half-space in wavenumber domain.
    %
    % Note that source strengh I and resistivity \rho only shift the
    % amplitude of the analytic solution such that dealing with unit
    % quantities is justified.
    %
    % The approach leads to k and w which approximate the bessel function
    % for arbitrary TX/RX offsets r.
    %
    % SYNTAX
    %   [k, w, n] = getInvFTBoerner(TX, RX, n)
    %
    % INPUT PARAMETER
    %   TX ... Matrix of TX position coordinates.
    %   RX ... Matrix of RX position coordinates.
    %
    % OUTPUT PARAMETER
    %   k ... Vector, containing quadrature nodes (spatial wavenumbers).
    %   w ... Vector, containing weights for quadrature summation.
    %   n ... Scalar, denoting number of k or w, respectively.
    %
    % REMARKS
    %   Also see Kemna A, Dissertation, 1999; LaBrecque et al 1996a
    
    
    %% Prepare parameter.
    
    % Get TX-RX offsets.
    r_TX2RX = [];
    for kk = 1:size (TX, 1)
        cur_TX2RX = sqrt((TX(kk,1) - RX(:,1)).^2 + (TX(kk,2) - RX(:,2)).^2);
        
        % Exclude TX == RX.
        cur_TX2RX(cur_TX2RX == 0) = [];
        
        % Summarize.
        r_TX2RX = [r_TX2RX; cur_TX2RX]; %#ok
    end
    r_TX2RX = unique(r_TX2RX);
    
    % Get critical wavenumber (transition of asymptotic behavior).
    k_crit = 1 / (2 * min(r_TX2RX));

    % Define number of wavenumbers for both parts of the asymptotic 
    % behavior of the analytic solution.
    % (Adapt unsymmetric distribution of nodes from Kemna)
    n_k_log = round(3*log(max(r_TX2RX) / min(r_TX2RX)));
    n_k_exp = 4;

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
    n = n_k_log + n_k_exp;
end

function [k, n] = getInvFTBing(TX, RX, n)
    % Calculate nodes and weights for quadrature approach by Bing Z. 1998.
    %
    % The values are based on the consideration of the bessel function
    % behaviour for variing wavenumbers w.r.t. the minimal and maximal
    % TX-RX offset.
    %
    % SYNTAX
    %   [k, n] = getInvFTBing(TX, RX)
    %
    % INPUT PARAMETER
    %   TX ... Matrix of TX position coordinates.
    %   RX ... Matrix of RX position coordinates.
    %   n  ... Scalar, denoting number of k or w, respectively.
    %
    % OUTPUT PARAMETER
    %   k   ... Vector, containing wavenumbers.
    %   n   ... Scalar, denoting number of k or w, respectively.
    
    %% Prepare parameter.
    
    % Get TX-RX offsets.
    r_TX2RX = [];
    for kk = 1:size (TX, 1)
        cur_TX2RX = sqrt((TX(kk,1) - RX(:,1)).^2 + (TX(kk,2) - RX(:,2)).^2);
        
        % Exclude TX == RX.
        cur_TX2RX(cur_TX2RX == 0) = [];
        
        % Summarize.
        r_TX2RX = [r_TX2RX; cur_TX2RX]; %#ok
    end
    r_TX2RX = unique(r_TX2RX);
    r_min = min(r_TX2RX);
    
    % Set lower bnd for f(k) and k.
    % (Derived from bessle function behavior for r_min = 0.5 and 
    % r_max = 500 -> we don't expect offsets beyond that range)
    u_min = 1e-4;
    k_min = 1e-3;
        
    %% Get initial k range.
    
    % Predefine a large range of wavenumbers.
    k_range = logspace(log10(k_min), 6, 301);

    % Define analytic solution for unit source over unit half-space.
    I = 1;
    rho = 2 * pi;
    u_HS = @(k, r) (I * rho) / (2 * pi) * besselk(0, k * r);

    % Get solutions spectra.
    u = arrayfun(@(k) u_HS(k, r_min), k_range);
    
    % Get k for u_max.
    k_max = k_range(find(u > u_min, 1, 'last'));
    
     % Get k from Boerner.
    k_crit = 1 / (2 * min(r_TX2RX));
    
    % Get range of k.
    k_2max = logspace(log10(k_crit), log10(k_max), 1 + round(n/2));
    k_2min = logspace(log10(k_crit), log10(k_min), n - round(n/2));
    k = [fliplr(k_2min), k_2max(2:end)].';
end

function [k, w, n] = getInvFTXu(TX, RX, n)
    % Calculate weights by Xu S. 2000 for quadrature nodes by Bing Z 1998.
    %
    % Approach uses bessel function K_0(k * r) behaviour for minimal and
    % maximal TX/RX offset to derive a range for k.
    % Use critical k from Boerner to derive a distribution for k.
    % Calculate the appropriate w by least-squares approach.
    %
    % The approach leads to optimized w which approximate the bessel
    % function best at a given set of TX/RX offsets r.
    %
    % SYNTAX
    %   [k, n] = getInvFTXu(TX, RX)
    %
    % INPUT PARAMETER
    %   TX ... Matrix of TX position coordinates.
    %   RX ... Matrix of RX position coordinates.
    %   n  ... Scalar, denoting number of k or w, respectively.
    %
    % OUTPUT PARAMETER
    %   k   ... Vector, containing wavenumbers.
    %   n   ... Scalar, denoting number of k or w, respectively.
    
    %% Prepare parameter.
    
    % Get TX-RX offsets.
    r_TX2RX = [];
    for kk = 1:size (TX, 1)
        cur_TX2RX = sqrt((TX(kk,1) - RX(:,1)).^2 + (TX(kk,2) - RX(:,2)).^2);
        
        % Exclude TX == RX.
        cur_TX2RX(cur_TX2RX == 0) = [];
        
        % Summarize.
        r_TX2RX = [r_TX2RX; cur_TX2RX]; %#ok
    end
    r_TX2RX = unique(r_TX2RX);
    r_min = min(r_TX2RX);
    
    % Set lower bnd for f(k) and k.
    % (Derived from bessle function behavior for r_min = 0.5 and 
    % r_max = 500 -> we don't expect offsets beyond that range)
    u_min = 1e-4;
    k_min = 1e-3;
        
    %% Get initial k range.
    
    % Predefine a large range of wavenumbers.
    k_range = logspace(log10(k_min), 6, 301);

    % Define analytic solution for unit source over unit half-space.
    I = 1;
    rho = 2 * pi;
    u_HS = @(k, r) (I * rho) / (2 * pi) * besselk(0, k * r);

    % Get solutions spectra.
    u = arrayfun(@(k) u_HS(k, r_min), k_range);
    
    % Get k for u_max.
    k_max = k_range(find(u > u_min, 1, 'last'));
    
     % Get k from Boerner.
    k_crit = 1 / (2 * min(r_TX2RX));
    
    % Get range of k.
    k_2max = logspace(log10(k_crit), log10(k_max), 1 + round(n/2));
    k_2min = logspace(log10(k_crit), log10(k_min), n - round(n/2));
    k = [fliplr(k_2min), k_2max(2:end)].';  
  
    %% Get corresponding w by linear least-squares.

    % The solution of point source over a homogeneous half-space 
    %   i) (\rho * I) / (2 * pi * abs(r))
    % is given by the Fourier transform in wavenumberdomain by
    %   ii) (2 / pi) \int_{0}^{\inf} (I * rho) / (2 * pi) * besselk(0, k * r)
    %
    % Hence, i) is approximated by quadrature approach applied on ii)
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
    %
    % Procedure: 
    %   1) take k -> Bing-approach
    %   2) optimize for w
    
    % Set up linear system for estimation of weights.
    A = getSys(u_HS, k, r_TX2RX);
    b = 1 + zeros(length(r_TX2RX), 1);
    
    % Obtain SVD.
    [U, S, V] = svd(A);
    sing_val = diag(S);
    pseudo_rank = find(sing_val > 1e-3, 1, 'last');
    invA = V(:,1:pseudo_rank) * diag(1./sing_val(1:pseudo_rank)) * U(:, 1:pseudo_rank).';

    % Solve minimization by Pseudoinverse.
    % (Add prefactor to derive weights,  equivalent to Boerner and Bing.)
    w = (pi / 2) * invA * b;

    %% Optimize k by nonlinear least-squares.
%{   
    % Iterate.
    iter_max = 100;
    res = zeros(iter_max, 1);
    for kk = 1:iter_max
        % Get Jacobian.
        A0 = getSys(u_HS, k, r_TX2RX);
        J = getSens(u_HS, k, r_TX2RX, A0, w);
        
        % Obtain SVD.
        [U, S, V] = svd(J);
        sing_val = diag(S);
        pseudo_rank = find(sing_val > sing_val(1) * 1e-3, 1, 'last');
        invJ = V(:,1:pseudo_rank) * diag(1./sing_val(1:pseudo_rank)) * U(:, 1:pseudo_rank).';
        
        % Set up rhs.
        b = ones(length(r_TX2RX), 1) - (A0 * w);
        res(kk) = norm(b);
        
        % Solve minimization by Pseudoinverse.
        % (Add prefactor to derive weights,  equivalent to Boerner and Bing.)
        dk = invJ * b;
        
        % Get model update.
        k = k + dk;
        
%         if norm(dk) < norm(k)*1e-6
%             break;
%         end
    end
%}
    %% Optimize w and k by nonlinear least-squares.
%{   
    % Nonlinear least-squares approach:
    %   ||I - F(k_j,w_j)||_2^2 -> min_[k_j,w_j]
    %
    % Linearization:
    %   F(k_j,w_j) ~= F(k0_j,w0_j) + [\d_F(k0_j,w0_j) / d_k_j; * [d_k_j,d_w_j]
    %                                 \d_F(k0_j,w0_j) / d_w_j]
    %
    % Gauss-Newton problem:
    %   ||I - F(k0_j,w0_j) - [\d_F(k0_j,w0_j) / d_k_j; * [d_k_j,d_w_j] ||_2^2 -> min_[d_k_j,d_w_j]
    %                         \d_F(k0_j,w0_j) / d_w_j] 
    %     |______ ______|  - |_________ _____________| * |_____ _____|
    %            V                     V                       V
    %            b         -     J(k0_j, w0_j)         *      d_m
    
    % Set up optimization for m = [k; w].
    w = 1 + 0*k;
    m = log([k; w]);
    
    iter = 1;
    iter_max = 300;
    res = zeros(iter_max, 1);
    while iter <= iter_max
        
        % Get derivative w.r.t. w at k_start: 
        % (== A since dependencies are linear!)
        S_w = getSys(u_HS, exp(m(1:n)), r_TX2RX);
        
        % Get derivative w.r.t. k: (using perturbation method)
        S_k = getSens(u_HS, exp(m(1:n)), r_TX2RX, S_w, exp(m(n + 1:end)));
        
        % Apply chain rule for darivatives w.r.t. log(k, w):
        S_w = S_w * diag(exp(m(n + 1:end)));
        S_k = S_k * diag(exp(m(1:n)));
        
        % Summarize Jacobian.
        J = [S_k, S_w];
        
        % Set up rhs.
        b = ones(length(r_TX2RX), 1) - (S_w * m(n + 1:end));

        % Get current residuum.
        res(iter) = norm(b);
        
%         % Obtain SVD.
%         [U, S, V] = svd(J);
%         sing_val = diag(S);
%         pseudo_rank = find(sing_val > sing_val(1)*1e-4, 1, 'last');
%         invJ = V(:,1:pseudo_rank) * diag(1./sing_val(1:pseudo_rank)) * U(:, 1:pseudo_rank).';
%         
%         % Solve.
%         d_m = invJ * b;
        
        % Solve.
        JTJ = J.'*J;
        d_m = (JTJ + 1e-6*eye(size(JTJ))) \ (J.'*b);
        
        % Get model update.
        m = m + d_m;
        iter = iter + 1;
    end
    
    % Exclude wavenumbers and weights.
    % (Add prefactor to derive weights, equivalent to Boerner and Bing.)
    k = exp(m(1:n));
    w = (pi / 2) * exp(m(n + 1:end));
%}
end

function A = getSys(fun, k, r)
    % Calculates system matrix of linear problem w.r.t. w.
   
    A = zeros(length(r), length(k));
    for ii = 1:length(r)
        A(ii,:) = r(ii) * arrayfun(@(x) fun(x,r(ii)), k);
    end
end

function S = getSens(fun, k, r, A0, w) %#ok
    % Calculates d_(F(k,w)) / d_k by perturbation.
    %
    % F(k,w) := r_l * A_lj(k_j, r_l) * w_j
    % -> getSys provides: r_l * A_lj(k_j, r_l)
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
        
        % Calculate colums of sensitivity.
        S(:, ii)= ((A_dist - A0) * w) / (k(ii).*noise);
    end
end

function [] = compare_w_k(TX, RX, n_k) %#ok
    % Get w and k:
    [k_xu, w_xu, ~] = getInvFTXu(TX, RX, n_k);
    [k_rub, w_rub, ~] = getInvFTBoerner(TX, RX);

    % Get TX-RX offsets.
    r_TX2RX = [];
    for kk = 1:size (TX, 1)
        cur_TX2RX = sqrt((TX(kk,1) - RX(:,1)).^2 + (TX(kk,2) - RX(:,2)).^2);

        % Exclude TX == RX.
        cur_TX2RX(cur_TX2RX == 0) = [];

        % Summarize.
        r_TX2RX = [r_TX2RX; cur_TX2RX]; %#ok
    end
    r_TX2RX = unique(r_TX2RX);
    n_r = length(r_TX2RX);

    % Calculate approximation.
    [r_approx_xu, r_approx_rub] = deal(zeros(n_r, 1));
    for kk = 1:n_r
        r_approx_xu(kk) = sum((2 / pi)*w_xu .* besselk(0, k_xu * r_TX2RX(kk)));
        r_approx_rub(kk) = sum((2 / pi)*w_rub .* besselk(0, k_rub * r_TX2RX(kk)));
    end

    % Compare asymptotic with approximation.
    figure(100)
    subplot(2, 1, 1)
        plot(1:n_r, 1./r_TX2RX, '-k', ...
             1:n_r, r_approx_xu, '-or', ...
             1:n_r, r_approx_rub, '-xb');
        legend('1/r_l', '1/r_{xu, l}', '1/r_{rub, l}');
        ylabel('f(r)');
        xlim([0.5, n_r+0.5]);
        xticklabels([]);
    subplot(2, 1, 2)
        plot(1:n_r, 100*(r_approx_xu - 1./r_TX2RX) / 1./r_TX2RX, '-ob', ...
             1:n_r, 100*(r_approx_rub - 1./r_TX2RX) / 1./r_TX2RX, '-vr');
        legend('1/r_{xu, l}', '1/r_{rub, l}');
        xlim([0.5, n_r+0.5]);
        xticks(1:n_r);
        xticklabels(r_TX2RX);
        xlabel('r');
        ylabel('rel. error');
end