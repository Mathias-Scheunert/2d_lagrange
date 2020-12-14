function u = unit(n, i, make_sparse)
    % Constructs i-th unit vector in n dimensions.
    %
    % SYNTAX
    %
    %   u = unit(n, i[, make_sparse])
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   n           ... Length of the resulting unit vector or number of
    %                   rows for the resulting matrix of unit vectors.
    %   i           ... Denotes which unit vector to construct, i.e. which
    %                   row of the otherwise all-zeros column vector should
    %                   be one. Constructs multiple unit vectors if 'i' is
    %                   a vector in which case each column of 'u'
    %                   corresponds to one entry in 'i'.
    %   make_sparse ... Logical scalar defaulting to false. If true, the
    %                   resulting vector or matrix 'u' will be sparse and
    %                   dense otherwise.
    %   u           ... Resulting unit column vector or matrix of unit
    %                   column vectors.

    if nargin < 3
        make_sparse = false;
    else
        assert(isscalar(make_sparse) && islogical(make_sparse), ...
            'Expected argument ''make_sparse'' to be a logical scalar.');
    end

    ni = length(i);
    u = accumarray({i, 1:ni}, 1, [n, ni], [], [], make_sparse);
end
