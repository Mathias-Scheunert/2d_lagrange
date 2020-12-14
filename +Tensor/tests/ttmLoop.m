function [res, used] = ttmLoop(ten, mat, dim)
    % TODO: Help.

    % Check arguments.
    assert(ndims(ten) <= 3, ...
        'Expected ''ten'' to be a tensor of order three.');
    assert(ismatrix(mat), ...
        'Expected ''mat'' to be a matrix.');
    assert(isscalar(dim) && any(dim == [1, 2, 3]), ...
        'Expected ''dim'' to be one of [1, 2, 3].');
    assert(size(ten, dim) == size(mat, 1), ...
        'Incompatible size between ''ten'' and ''mat''.');

    % Get sizes.
    [m, n, p] = size(ten);
    [~, q] = size(mat);

    % Compute.
    if dim == 1
        used = false(m, q, n, p);
        res = zeros(q, n, p);
        for k = 1:p
            for j = 1:n
                for i = 1:q
                    u = asRow(ten(:, j, k));
                    v = mat(:, i);
                    res(i, j, k) = u * v;
                    used(:, i, j, k) = u.' & v;
                end
            end
        end
    elseif dim == 2
        used = false(m, n, q, p);
        res = zeros(m, q, p);
        for k = 1:p
            for j = 1:q
                for i = 1:m
                    u = asRow(ten(i, :, k));
                    v = mat(:, j);
                    res(i, j, k) = u * v;
                    used(i, :, j, k) = u.' & v;
                end
            end
        end
    else % dim == 3
        used = false(m, n, p, q);
        res = zeros(m, n, q);
        for k = 1:q
            for j = 1:n
                for i = 1:m
                    u = asRow(ten(i, j, :));
                    v = mat(:, k);
                    res(i, j, k) = u * v;
                    used(i, j, :, k) = u.' & v;
                end
            end
        end
    end
end
