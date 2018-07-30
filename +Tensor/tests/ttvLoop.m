function res = ttvLoop(ten, vec, dim)
    % TODO: Help.
    
    % Check arguments.
    assert(ndims(ten) <= 3, ...
        'Expected ''ten'' to be a tensor of order three.');
    assert(iscolumn(vec), ...
        'Expected ''vec'' to be a column vector.');
    assert(isscalar(dim) && any(dim == [1, 2, 3]), ...
        'Expected ''dim'' to be one of [1, 2, 3].');
    assert(size(ten, dim) == length(vec), ...
        'Incompatible size between ''ten'' and ''vec''.');
    
    % Get sizes.
    [m, n, p] = size(ten);
    
    % Compute.
    if dim == 1
        res = zeros(n, p);
        for j = 1:p
            for i = 1:n
                res(i, j) = asRow(ten(:, i, j)) * vec;
            end
        end
    elseif dim == 2
        res = zeros(m, p);
        for j = 1:p
            for i = 1:m
                res(i, j) = asRow(ten(i, :, j)) * vec;
            end
        end
    else % dim == 3
        res = zeros(m, n);
        for j = 1:n
            for i = 1:m
                res(i, j) = asRow(ten(i, j, :)) * vec;
            end
        end
    end
end
