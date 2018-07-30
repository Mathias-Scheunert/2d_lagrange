% Tests multiplication of sparse tensor classes.

% Tensor class.
ten_class = @Tensor.Tensor3Coord;

% Specify size.
n = 3;
m = 5;

% Construct index vectors.
i = mod((0:m - 1), n) + 1;
k = 0 + 1 * (1:m);
iik = [i(:), i(:), k(:)];

% Genereate pseudo-random values.
v_rs = RandStream('mt19937ar', 'Seed', 20120515);
v = rand(v_rs, size(iik, 1), 1);

% Construct sparse tensor.
t = ten_class([n, n, m], iik, v);
t_full = full(t);
display(t_full);

% Test tensor-vector product.
fprintf('\nTesting ''ttv'':\n');
v3 = pick(1, unit(m, 1, true), ones(m, 1));
tv3_coor = full(ttv(t, v3, 3));
tv3_loop = ttvLoop(t_full, full(v3), 3);
display(tv3_coor);
display(tv3_loop);

% Test tensor-matrix product.
fprintf('\nTesting ''ttm'':\n');
m3 = sparse((1:m).', [1, 2, 1, 2, 1].', 1, 5, 2);
tm3_coor = full(ttm(t, m3, 3));
tm3_loop = ttmLoop(t_full, full(m3), 3);
display(tm3_coor);
display(tm3_loop);

% Verify equality.
assert(isequal(tv3_coor, tv3_loop), ...
    'Mismatch between ''tv3_coor'' and ''tv3_loop''.');
assert(isequal(tm3_coor, tm3_loop), ...
    'Mismatch between ''tm3_coor'' and ''tm3_loop''.');
