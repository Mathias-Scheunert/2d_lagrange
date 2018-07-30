function [res_i, res_v, used] = ttmCoordDebug(ten_i, ten_m, ten_v, mat_i, mat_j, mat_m, mat_v, dim, sz, mat_sz)
    % TODO: Help.
    
    % Debug version ('ttmCoordCount' & 'ttmCoordApply') to be inserted in
    % 'Tensor3Coord.ttm'.
    %{
    [res_i, res_v, used] = ttmCoordDebug(...
        ten_i, ten_m, ten_v, mat_i, mat_j, mat_m, mat_v, dim, sz, mat_sz);
    %}
    
    used_sz = [sz(1:dim - 1), mat_sz, sz(dim + 1:3)];
    used = false(used_sz);
    res_i = zeros(0, 3);
    res_v = zeros(0, 1);
    num_nnz = 0;
    ten_num = size(ten_i, 1);
    ten_cur = 1;
    ten_off = ten_cur;
    mat_num = size(mat_i, 1);
    idx_match = false;
    %
    while ten_cur <= ten_num
        mat_cur = 1;
        mat_off = mat_cur;
        %
        while mat_cur <= mat_num
            if ten_i(ten_cur, dim) == mat_i(mat_cur)
                if idx_match
                    old = res_v(num_nnz, 1);
                else
                    idx_match = true;
                    num_nnz = num_nnz + 1;
                    old = 0;
                    res_i(num_nnz, :) = ten_i(ten_cur, :);
                    res_i(num_nnz, dim) = mat_j(mat_cur);
                end
                res_v(num_nnz, 1) = old + ...
                    ten_v(ten_cur) .* mat_v(mat_cur);
                %
                used_idx = num2cell([...
                    ten_i(ten_cur, 1:dim), ...
                    mat_j(mat_cur), ...
                    ten_i(ten_cur, dim + 1:3)]);
                used(used_idx{:}) = true;
                %
                ten_cur = ten_cur + 1;
                mat_cur = mat_cur + 1;
            elseif ten_i(ten_cur, dim) < mat_i(mat_cur)
                ten_cur = ten_cur + 1;
            else % ten_i(ten_cur, dim) > mat_i(mat_cur)
                mat_cur = mat_cur + 1;
            end
            %
            if mat_m(mat_cur) && mat_cur ~= mat_off
                idx_match = false;
                mat_off = mat_cur;
                ten_cur = ten_off;
            end
            %
            if ten_m(ten_cur) && ten_cur ~= ten_off
                idx_match = false;
                while ~mat_m(mat_cur) || mat_cur == mat_off
                    mat_cur = mat_cur + 1;
                end
                mat_off = mat_cur;
                ten_cur = ten_off;
            end
        end
        %
        idx_match = false;
        while ~ten_m(ten_cur) || ten_cur == ten_off
            ten_cur = ten_cur + 1;
        end
        ten_off = ten_cur;
    end
end
