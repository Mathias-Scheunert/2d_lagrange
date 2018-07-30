% Tests size-related queries for sparse tensor classes.

% Tensor class.
ten_class = @Tensor.Tensor3Coord;

% Create all-zero tensor.
t = ten_class([5, 5, 2])

% Query various size-related functions.
ie = isempty(t)
nd = ndims(t)
ne = numel(t)
sz = size(t)
d1 = size(t, 1)
d2 = size(t, 2)
d3 = size(t, 3)
d4 = size(t, 4)
[s2{1:2}] = size(t)
[s3{1:3}] = size(t)
[s4{1:4}] = size(t)

%#ok<*SNASGU,*NOPTS>
