function [fit, last_rate] = getConvRates(errtype, num_dofs, errs, verbosity)
    % Calculates convergence rates.
    %
    % SYNTAX
    %   [fit, last_rate] = getConvRates(errtype, num_dofs, errs[, verbosity])
    %
    % INPUT PARAMETERS
    %   errtype  ... String, denoting the error.
    %   num_dofs ... Cell [n x 1], containing the DOF for all calculated n
    %                orders of FE.
    %   errs     ... Cell [n x 1], containing the respective errors for all
    %                 calculated n orders of FE.
    %
    % OUTPUT PARAMETERS
    %   fit       ... Scalar, denoting the mean slope of the error curve
    %                 w.r.t. the DOF.
    %   last_rate ... Scalar, denoting the slope of the last two entries of
    %                 the error curve w.r.t. the DOF.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % COPYRIGHT
    %   Originally written by Jan Blechta (CurlCurl-Toolbox).
    
    %% Check imput.
    
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    % Set problem dimension.
    dim = 2;
    
    %% Calculate and print.
    
    if verbosity
        fprintf('%s-error vs (#dofs)^(1/%d): LS fit, last rate\n', ...
            errtype, dim);
    end
    
    [fit, last_rate] = deal(cell(size(num_dofs, 1), 1));
    for order = 1:size(num_dofs, 1)
        % Calculate mean convergance rate.
        cur_fit = polyfit(log10(num_dofs{order}), dim*log10(errs{order}), 1);
        fit{order} = cur_fit(1);
        
        % Calculate rate of two last refninement steps.
        last_rate{order} = dim * log(errs{order}(end)/errs{order}(end-1)) ...
          / log(num_dofs{order}(end)/num_dofs{order}(end-1));
      
        if verbosity
            fprintf('    order %d, LS %1.2f, last %1.2f\n', ...
                order, fit{order}, last_rate{order});
        end
    end
end