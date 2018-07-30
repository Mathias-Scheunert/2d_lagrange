% Tests multiplication of sparse tensor class Tensor3Coord thoroughly.

% Sizes of tensor and matrix as well as density.
m = 29;
n = 23;
p = 19;
q = 17;
T_density = 0.1;
M_density = 0.1;
val_choice = 1;

% Initialize random number generator.
rs = RandStream('mt19937ar', 'Seed', 20120522);

% Derive sizes and number of nonzeros.
T_sz = [m, n, p];
M_sz = arrayfun(@(x) {[x, q]}, T_sz);
T_nnz = floor(T_density * prod(T_sz));
M_nnz = cellfun(@(X_sz) floor(M_density * prod(X_sz)), M_sz);

% Create tensor and matrices.
T_index = arrayfun(@(x) {rs.randi(x, T_nnz, 1)}, T_sz);
T_value = pick(val_choice, rs.rand(T_nnz, 1), ones(T_nnz, 1));
T = Tensor3Coord(T_sz, cell2mat(T_index), T_value);
T_full = full(T);
M = cell(1, 3);
for i = 1:3
    M_index = arrayfun(@(x) {rs.randi(x, M_nnz(i), 1)}, M_sz{i});
    M_value = pick(val_choice, rs.rand(M_nnz(i), 1), ones(M_nnz(i), 1));
    M{i} = sparse(M_index{:}, M_value, M_sz{i}(1), M_sz{i}(2));
end
M_full = cellfun(@(X) {full(X)}, M);

% Test.
for i = 1:3
    fprintf('[dim = %d]\n', i);
    %{
    [B_coor, used_coor] = ttm(T, M{i}, i);
    [B_loop, used_loop] = ttmLoop(T_full, M_full{i}, i);
    all_used = [nnz(used_coor), nnz(used_loop)];
    fprintf('  all_used = [%8d, %8d]\n', all_used);
    %}
    B_coor = ttm(T, M{i}, i);
    B_full = full(B_coor);
    B_loop = ttmLoop(T_full, M_full{i}, i);
    all_nnz = [nnz(B_full), nnz(B_loop)];
    fprintf('  all_nnz  = [%8d, %8d]\n', all_nnz);
    fprintf('  equal    =  %8d\n', isequal(B_full, B_loop));
    fprintf('  err_norm =  %.6e\n', norm(B_full(:) - B_loop(:)));
end
