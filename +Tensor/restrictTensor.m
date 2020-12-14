function X = restrictTensor(X, Null)
    % Restricts a tensor to a subset of its rows and columns.
    %
    % SYNTAX
    %
    %   X = restrictTensor(X, Null)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   X    ... On input, a sparse tensor of size n x n x k. On output, a
    %            sparse tensor of size m x m x k and of the same class as
    %            the input.
    %   Null ... Sparse matrix of size n x m, logical vector of length n
    %            with m entries being true, or a MatrixRestriction object.
    %            If provided as a sparse matrix, all entries are expected
    %            to be either zero or one and no row or column is allowed
    %            to have more than one nonzero element.
    %
    % REMARKS
    %
    %   The tensor 'X' is restricted along the first and second dimension
    %   according to the given restriction 'Null' while the third dimension
    %   remains unaltered. The effect is as if 'restrictMatrix' was applied
    %   to every slice 'X(:, :, i)' separately and the resulting smaller
    %   slices again concatenated along the third dimension. Another
    %   (albeit quite inefficient) way to express this would be:
    %
    %     X = ttm(ttm(X, Null, 1), Null, 2)
    %
    % See also restrictMatrix.

    % Verify that 'X' is indeed a tensor.
%     assert(isa(X, 'Tensor3Base') | isa(X, 'Tensor3Coord'), ...
%         'sputil:restrictTensor:NotTensorClass', [...
%         'Expected ''X'' to be an instance of ''Tensor3Base'' or one ', ...
%         'of its subclasses.']);

    % Get the matrix restriction and its representation as a mask.
    Null = Tensor.MatrixRestriction(Null);
    %Null = MatrixRestriction(Null);
    mask = getMask(Null);

    % Get size of matrix restriction and tensor. Check for compatibility.
    [num_row, num_col] = size(Null);
    sz = size(X);
    assert(all(sz(1:2) == num_row), ...
        'sputil:restrictTensor:IncompatibleSize', [...
        'Expected the size of the first two dimensions of ''X'' to ', ...
        'match the number of rows in ''Null''.']);

    % Disassemble tensor.
    class_fn = str2func(class(X));
    [index, value] = find(X);

    % Eliminate indices and values not in mask.
    valid = mask(index(:, 1)) & mask(index(:, 2));
    index = index(valid, :);
    value = value(valid);

    % Renumber first two dimensions and adjust size.
    map = zeros(num_row, 1);
    map(mask) = 1:num_col;
    index(:, 1:2) = map(index(:, 1:2));
    sz(1:2) = num_col;

    % Reassemble restricted tensor.
    X = class_fn(sz, index, value);
end
