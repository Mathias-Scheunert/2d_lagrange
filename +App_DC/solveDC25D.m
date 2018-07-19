function u = solveDC25D(fe, sol, FT, verbosity)
    % Apply quadrature strategies for inverse spatial Fourirer trafo.
    %
    % SYNTAX
    %   u = solveDC25D(fe, sol, FT, verbosity)
    %
    % INPUT PARAMETER
    %   fe  ... Struct, including all information to set up Lagrange 
    %           FE, as well as the linear system components.
    %   sol ... Cell [n x 1] of structs for the n wavenumbers.
    %           Each containing the information of the DC problem to be 
    %           solved numerically, i.e. rhs vector, system matrix, 
    %           interpolation operator.
    %   FT  ... Struct, containing the weights and or nodes and the
    %           quadrature type for inverse spatial FT.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   u ... Solution vector for 2.5D DC Lagrange FE problem.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(iscell(sol), ...
        'sol - Cell, containing sol structs for each wavenumber, expected.');
    assert(isstruct(FT) && all(isfield(FT, {'n', 'type'})), ...
        'FT - struct, including all information for inverse spatial FT, expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    %% Solve fwd problems.

    if verbosity
        fprintf('Solve linear systems for all wavenumbers ... '); 
    end
    u_2D = cell(FT.n, 1);
    for jj = 1:(FT.n)
        u_2D{jj} = Fe.solveFwd(sol{jj}, fe);
    end
    if verbosity
        fprintf('done.\n'); 
    end
    
    %% Apply numerical integration over wavenumber domain.
    
    if verbosity
        fprintf('Apply integration over wavenumber domain ... '); 
    end
    switch FT.type
        case 'Boerner'
            % Add quadrature weights.
            u = cellfun(@(x, y) {x * y}, u_2D, num2cell(FT.w));

            % Sum up solutions (apply quadrature).
            u = (2 / pi) .* sum(cat(3, u{:}), 3);
            
        case 'Bing'
            % Sum up solutions (apply quadrature).
            u_rect = (2  * FT.k(1) / pi) * u_2D{1};
            u_log = 0 * u_rect;
            for ii = 2:(length(FT.k) - 2)
                u_log = u_log + (...
                    (2 / pi) * (FT.k(ii + 1) - FT.k(ii)) * ...
                    (u_2D{ii} - u_2D{ii + 1}) ./ ...
                    (log(u_2D{ii} ./ u_2D{ii + 1}))...
                        );
            end
            u_exp = ((2 * u_2D{end}) * (FT.k(end) - FT.k(end - 1))) ./ ...
                (pi * log(u_2D{end - 1} ./ u_2D{end}));
            u = u_rect + u_log + u_exp;
            
        case 'Xu'
            % Add quadrature weights.
            u = cellfun(@(x, y) {x * y}, u_2D, num2cell(FT.w));

            % Sum up solutions (apply quadrature).
            u = sum(cat(3, u{:}), 3);
        otherwise
            error('Unknown type.');
    end
    if verbosity
        fprintf('done.\n');
        fprintf('... DC 2.5D problem solved.\n \n');
    end
end