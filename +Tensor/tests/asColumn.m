function x = asColumn(x)
    % Reshapes an array into a column vector.
    %
    % SYNTAX
    %
    %   x = asColumn(x)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   x ... On input, an arbitrary MATLAB array. On output, the same
    %         array reshaped into a column vector.
    %
    % REMARKS
    %
    %   The effect of this function is identical to the expression:
    %
    %     x = x(:)
    %
    %   However, it can be used in contexts where the array is not directly
    %   accessed by name and thus the indexing notation is not available.
    %   Possible usage areas are the function parameter of 'cellfun' or
    %   direct reshaping of a function's return value.
    %
    % See also asRow.

    x = reshape(x, numel(x), 1);
end
