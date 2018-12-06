function [] = createGmsh_Ralph(mn1, mn2, ab_half)
    % Construct test mesh for gmsh.
    
    % Set gmsh parameter.
    exp_fac = pick(1, 5, 10);          % extension of domain w.r.t. AB/2
    % -> increase to improve solution for large AB/2!
    h_at_TXRX = pick(1, 0.5);          % element size at TX/RX points
    
    % Set geopmetry.
    d = 5;                             % width of the vertical dike    
    offset = -0.1;                     % horz shift of configuration midpoint
    
    %% Define points describing geometry.

    % Domain (considering clock-wise sorting).
    bnd_x = max(ab_half) * exp_fac;    % define outermost model boundary in x
    bnd_z = bnd_x;                     % define outermost model boundary in z
    h_domain = bnd_x/10;

    p_domain = [-bnd_x, 0, 0, h_domain;
                 bnd_x, 0, 0, h_domain;
                 bnd_x, 0, bnd_z, h_domain;
                -bnd_x, 0, bnd_z, h_domain]; 

    % Dike (considering clock-wise sorting).
    h_dike_ymin = h_at_TXRX;
    h_dike_ymax = h_domain/2;
    horz_shift = pick(1, 0, bnd_x, -bnd_x/2); % turns vertical dike to dipping dike
                                              % Note: geometry is NOT preserved

    p_dike = [0, 0, 0, h_dike_ymin;
              d, 0, 0, h_dike_ymin;
              horz_shift + d, 0, bnd_z, h_dike_ymax;
              horz_shift + 0, 0, bnd_z, h_dike_ymax];

    % Electrodes MN
    h_mn = h_at_TXRX;

    p_mn = [offset - (mn1/2), 0, 0, h_mn;
            offset + (mn1/2), 0, 0, h_mn;
            offset - (mn2/2), 0, 0, h_mn;
            offset + (mn2/2), 0, 0, h_mn];

    % Electrodes AB
    n_ab = length(ab_half);
    h_ab = h_mn;

    p_ab = zeros(2*n_ab, 4);
    p_a_x = (offset - ab_half).';
    p_b_x = (offset + ab_half).';
    p_ab(:,1) = [p_a_x; p_b_x];
    p_ab(:,2) = 0;
    p_ab(:,3) = 0;
    p_ab(:,4) = h_ab;

    %% Summarize and sort to a global point list.

    % Combine all surface points and sort w.r.t. x-coordinate.
    p_surface = [p_domain(1:2,:);
                 p_dike(1:2,:);
                 p_mn;
                 p_ab];
    p_surface = sortrows(p_surface, 1);

    % Add points in subsurface (considering clock-wise sorting).
    p_subsurface = [p_domain(3,:);
                    p_dike(3:end,:);
                    p_domain(4,:)];

    % Summarize points.
    p_all = [p_surface;
             p_subsurface];
    n_p_all = size(p_all, 1);

    %% Get global indices for specific points.

    % Dike.
    [~, idx_dik_l_top] = ismember(p_all, p_dike(1,:), 'rows');
    idx_dik_l_top = find(idx_dik_l_top);
    [~, idx_dik_r_top] = ismember(p_all, p_dike(2,:), 'rows');
    idx_dik_r_top = find(idx_dik_r_top);
    [~, idx_dik_r_bot] = ismember(p_all, p_dike(3,:), 'rows');
    idx_dik_r_bot = find(idx_dik_r_bot);
    [~, idx_dik_l_bot] = ismember(p_all, p_dike(4,:), 'rows');
    idx_dik_l_bot = find(idx_dik_l_bot);

    % Domain.
    [~, idx_dom_l_top] = ismember(p_all, p_domain(1,:), 'rows');
    idx_dom_l_top = find(idx_dom_l_top);
    [~, idx_dom_r_top] = ismember(p_all, p_domain(2,:), 'rows');
    idx_dom_r_top = find(idx_dom_r_top);
    [~, idx_dom_r_bot] = ismember(p_all, p_domain(3,:), 'rows');
    idx_dom_r_bot = find(idx_dom_r_bot);
    [~, idx_dom_l_bot] = ismember(p_all, p_domain(4,:), 'rows');
    idx_dom_l_bot = find(idx_dom_l_bot);

    %% Define gmsh straight line entities.

    % Set list of point indices describing straight lines at y min/max.
    pt_dom1_ymin = idx_dom_l_top:idx_dik_l_top;
    pt_dom2_ymin = idx_dik_l_top:idx_dik_r_top;
    pt_dom3_ymin = idx_dik_r_top:idx_dom_r_top;

    pt_dom1_ymax = idx_dik_l_bot:idx_dom_l_bot;
    pt_dom2_ymax = idx_dik_r_bot:idx_dik_l_bot;
    pt_dom3_ymax = idx_dom_r_bot:idx_dik_r_bot;

    % Set list of point indices describing vertical straight lines.
    pt_xmin = [idx_dom_l_top, idx_dom_l_bot];
    pt_xmax = [idx_dom_r_top, idx_dom_r_bot];

    pt_dik_xmin = [idx_dik_l_top, idx_dik_l_bot];
    pt_dik_xmax = [idx_dik_r_top, idx_dik_r_bot];

    % Set list of lines at y min/max.
    l_dom1_ymin = [pt_dom1_ymin(1:end-1).', pt_dom1_ymin(2:end).'];
    l_dom2_ymin = [pt_dom2_ymin(1:end-1).', pt_dom2_ymin(2:end).'];
    l_dom3_ymin = [pt_dom3_ymin(1:end-1).', pt_dom3_ymin(2:end).'];

    l_dom1_ymax = [pt_dom1_ymax(1:end-1).', pt_dom1_ymax(2:end).'];
    l_dom2_ymax = [pt_dom2_ymax(1:end-1).', pt_dom2_ymax(2:end).'];
    l_dom3_ymax = [pt_dom3_ymax(1:end-1).', pt_dom3_ymax(2:end).'];

    % Set list of vertical lines.
    l_xmin = [pt_xmin(1:end-1).', pt_xmin(2:end).'];
    l_xmax = [pt_xmax(1:end-1).', pt_xmax(2:end).'];

    l_dik_xmin = [pt_dik_xmin(1:end-1).', pt_dik_xmin(2:end).'];
    l_dik_xmax = [pt_dik_xmax(1:end-1).', pt_dik_xmax(2:end).'];

    %% Export all points.

    file_ID = fopen('test.geo', 'w');

    % Export geometrical point entities.
    fprintf(file_ID, '// Point list\n');
    for ii = 1:n_p_all
        fprintf(file_ID, 'Point(%d) = {%d, %d, %d, %d};\n', ii, p_all(ii,:));
    end
    % Get MN index in point list.
    [~, m1_in_p_all] = ismember(p_mn(1,:), p_all, 'rows');
    [~, n1_in_p_all] = ismember(p_mn(2,:), p_all, 'rows');
    [~, m2_in_p_all] = ismember(p_mn(3,:), p_all, 'rows');
    [~, n2_in_p_all] = ismember(p_mn(4,:), p_all, 'rows');
    % Export MN physical points.
    fprintf(file_ID, 'Physical Point("M1") = {%d};\n', m1_in_p_all);
    fprintf(file_ID, 'Physical Point("N1") = {%d};\n', n1_in_p_all);
    fprintf(file_ID, 'Physical Point("M2") = {%d};\n', m2_in_p_all);
    fprintf(file_ID, 'Physical Point("N2") = {%d};\n', n2_in_p_all);
    % Export AB physical points.
    for jj = 1:length(p_a_x)
        [~, cur_a_in_p_all] = ismember([p_a_x(jj), 0, 0, h_ab], p_all, 'rows');
        fprintf(file_ID, 'Physical Point("A%d") = {%d};\n', jj, cur_a_in_p_all);
        [~, cur_b_in_p_all] = ismember([p_b_x(jj), 0, 0, h_ab], p_all, 'rows');
        fprintf(file_ID, 'Physical Point("B%d") = {%d};\n', jj, cur_b_in_p_all);
    end

    %% Export all lines.

    fprintf(file_ID, '\n// Line list\n');

    % Export geometrical line entities at ymin.
    % -------->
    for ii = 1:size(l_dom1_ymin, 1)
        idx_l_dom1_ymin = size(l_dom1_ymin, 1);
        fprintf(file_ID, 'Line(%d) = {%d, %d};\n', ii, ...
            l_dom1_ymin(ii,1), l_dom1_ymin(ii,2));
    end
    loop_idx_dom1_ymin = 1:idx_l_dom1_ymin;

    for ii = 1:size(l_dom2_ymin, 1)
        idx_l_dom2_ymin = idx_l_dom1_ymin + size(l_dom2_ymin, 1);
        fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dom1_ymin + ii, ...
            l_dom2_ymin(ii,1), l_dom2_ymin(ii,2));
    end
    loop_idx_dom2_ymin = (idx_l_dom1_ymin + 1):idx_l_dom2_ymin;

    for ii = 1:size(l_dom3_ymin, 1)
        idx_l_dom3_ymin = idx_l_dom2_ymin + size(l_dom3_ymin, 1);
        fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dom2_ymin + ii, ...
            l_dom3_ymin(ii,1), l_dom3_ymin(ii,2));
    end
    loop_idx_dom3_ymin = (idx_l_dom2_ymin + 1):idx_l_dom3_ymin;

    % Export physical lines at ymin.
    fprintf(file_ID, ['Physical Line("ymin") = {', ...
        repmat('%d, ', 1, idx_l_dom3_ymin - 1), '%d};\n'], ...
        1:idx_l_dom3_ymin);
    % <--------

    % Export geometrical line entities at ymax.
    % -------->
    for ii = 1:size(l_dom1_ymax, 1)
        idx_l_dom1_ymax = size(l_dom1_ymax, 1);
        fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dom3_ymin + ii, ...
            l_dom1_ymax(ii,1), l_dom1_ymax(ii,2));
    end
    loop_idx_dom1_ymax = idx_l_dom3_ymin + (1:idx_l_dom1_ymax);

    for ii = 1:size(l_dom2_ymax, 1)
        idx_l_dom2_ymax = idx_l_dom1_ymax + size(l_dom2_ymax, 1);
        fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dom3_ymin + idx_l_dom1_ymax + ii, ...
            l_dom2_ymax(ii,1), l_dom2_ymax(ii,2));
    end
    loop_idx_dom2_ymax = idx_l_dom3_ymin + ((idx_l_dom1_ymax + 1):idx_l_dom2_ymax);

    for ii = 1:size(l_dom3_ymax, 1)
        idx_l_dom3_ymax = idx_l_dom2_ymax + size(l_dom3_ymax, 1);
        fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dom3_ymin + idx_l_dom2_ymax + ii, ...
            l_dom3_ymax(ii,1), l_dom3_ymax(ii,2));
    end
    loop_idx_dom3_ymax = idx_l_dom3_ymin + ((idx_l_dom2_ymax + 1):idx_l_dom3_ymax);

    % Export physical lines at ymax.
    fprintf(file_ID, ['Physical Line("ymax") = {', ...
        repmat('%d, ', 1, idx_l_dom3_ymax - 1), '%d};\n'], ...
        idx_l_dom3_ymin + (1:idx_l_dom3_ymax));
    % <--------

    % Export geometrical/physical line entities at xmin.
    idx_l_xmin = idx_l_dom3_ymin + idx_l_dom3_ymax + 1;
    fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_xmin, ...
            l_xmin(1), l_xmin(2));
    fprintf(file_ID, 'Physical Line("xmin") = {%d};\n', ...
        idx_l_dom3_ymin + idx_l_dom3_ymax + 1);

    % Export geometrical/physical line entities at xmax.
    idx_l_xmax = idx_l_dom3_ymin + idx_l_dom3_ymax + 2;
    fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_xmax, ...
            l_xmax(1), l_xmax(2));
    fprintf(file_ID, 'Physical Line("xmax") = {%d};\n', ...
        idx_l_dom3_ymin + idx_l_dom3_ymax + 2);

    % Export geometrical line/physical entities at dike x min/max.
    idx_l_dike_xmin = idx_l_dom3_ymin + idx_l_dom3_ymax + 3;
    fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dike_xmin, ...
            l_dik_xmin(1), l_dik_xmin(2));
    fprintf(file_ID, 'Physical Line("dike_xmin") = {%d};\n', ...
        idx_l_dom3_ymin + idx_l_dom3_ymax + 3);
    idx_l_dike_xmax = idx_l_dom3_ymin + idx_l_dom3_ymax + 4;
    fprintf(file_ID, 'Line(%d) = {%d, %d};\n', idx_l_dike_xmax, ...
            l_dik_xmax(1), l_dik_xmax(2));
    fprintf(file_ID, 'Physical Line("dike_xmax") = {%d};\n', ...
        idx_l_dom3_ymin + idx_l_dom3_ymax + 4);

    %% Export all surfaces.

    fprintf(file_ID, '\n// Surface list\n');
    % Export line loop for surface 1 = domain 1 (left of dike).
    l_loop_dom1 = sort([loop_idx_dom1_ymin, idx_l_dike_xmin, ...
                        loop_idx_dom1_ymax, -idx_l_xmin]);
    idx_loop_dom1 = 1;
    fprintf(file_ID, ['Line Loop(%d) = {', ...
        repmat('%d,', 1, length(l_loop_dom1) - 1), '%d};\n'], idx_loop_dom1, ...
        l_loop_dom1);
    fprintf(file_ID, 'Plane Surface(%d) = {%d};\n', idx_loop_dom1, idx_loop_dom1);
    fprintf(file_ID, 'Physical Surface("dom1") = {%d};\n', idx_loop_dom1);

    % Export line loop for surface 2 (dike).
    l_loop_dom2 = sort([loop_idx_dom2_ymin, idx_l_dike_xmax, ...
                        loop_idx_dom2_ymax, -idx_l_dike_xmin]);
    idx_loop_dom2 = 2;
    fprintf(file_ID, ['Line Loop(%d) = {', ...
        repmat('%d,', 1, length(l_loop_dom2) - 1), '%d};\n'], idx_loop_dom2, ...
        l_loop_dom2);
    fprintf(file_ID, 'Plane Surface(%d) = {%d};\n', idx_loop_dom2, idx_loop_dom2);
    fprintf(file_ID, 'Physical Surface("dom2") = {%d};\n', idx_loop_dom2);

    % Export line loop for surface 2 (right of dike).
    l_loop_dom3 = sort([loop_idx_dom3_ymin, idx_l_xmax, ...
                        loop_idx_dom3_ymax, -idx_l_dike_xmax]);
    idx_loop_dom3 = 3;
    fprintf(file_ID, ['Line Loop(%d) = {', ...
        repmat('%d,', 1, length(l_loop_dom3) - 1), '%d};\n'], idx_loop_dom3, ...
        l_loop_dom3);
    fprintf(file_ID, 'Plane Surface(%d) = {%d};\n', idx_loop_dom3, idx_loop_dom3);
    fprintf(file_ID, 'Physical Surface("dom3") = {%d};\n', idx_loop_dom3);

    fclose(file_ID);
end