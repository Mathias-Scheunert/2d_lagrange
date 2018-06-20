function const_fun = getConstFunction(val)
    % Get information for constant "function".
    %
    % This function mainly aims to mimic the structure of the other, more
    % complex, predefined functions to the trivial case of a constant.
    %
    % SYNTAX
    %
    %   const_fun = getConstFunction(val)
    %
    % INPUT PARAMETER
    %   val ... Scalar, denoting the amplitude.
    %
    % OUTPUT PARAMETER
    %
    %   const_fun ... Struct, containing function handle.
    
    %% Check input
    
    assert(isscalar(val), ...
        'val - scalar, denoting the amplitude, expected.');
    
    %% Define const-fuction.

    const_fun.f = @(x, y) val;
end