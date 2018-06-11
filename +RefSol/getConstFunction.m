function Const = getConstFunction(val)
    % Get information for constant function.
    %
    % SYNTAX
    %
    %   Const = getConstFunction(val)
    %
    % INPUT PARAMETER
    %   val ... Scalar, denoting the amplitude.
    %
    % OUTPUT PARAMETER
    %
    %   Const ... Struct, containing function handle.
    
    %% Check input
    
    assert(isscalar(val), ...
        'val - scalar, denoting the amplitude, expected.');
    
    %% Define const-fuction.

    Const.f = @(x, y) val;
end