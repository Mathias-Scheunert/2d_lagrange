function mesh = createChecker(anz,sizes,group,start_x,point_ele)
% creates mesh with checkered interior,
% assumes electrodes range from -60 to 60 if not specified otherwise
% with a distance of 1m
% anz: m x n amount of checker fields
% size: size of checker fields (uniform or varied with x and depth), e.g for 3
% x 4 checkers: 7 values , first values: y, uniform sizes when only 2
% values
%start_x: where x value where checkers begin, if empty checker field is
%centered
%      __________________________________________
%     |           |__0__|___1___|                 |
%     |           |__1__|___0___|                 |
%     |           |__0__|___1___|                 |
%     |                                           |
%     |                  0                        |
%     |___________________________________________|
%
 %% Check Gmsh version.

    gmsh_path = dir('**/gmsh.exe');
    if isempty(gmsh_path)
        error('Gmsh executable could not be found in path.');
    end
%% Create Gmsh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gmsh_file = 'checker';
%    input_info = createGmshInput(anz,sizes,group,start_x);
%    writeGEO([gmsh_file, '.geo'], input_info)
    %  Note that the 2D grid will be referred to gmesh x and z coordinate:
    %  0,0 ------>
    %      |      x = gmesh-x
    %      |
    %      |
    %      v
    %       y = gmesh-z
    min_x   = -10000;
    max_x   = 10000;
    max_y   = 10000;
    bnd = [min_x, max_x, 0, max_y];

%% Electrode points


    if isempty(point_ele)
        length_ele = 121;
        point_ele = [(-60:60)',zeros(121,1)];
    else


        length_ele = length(point_ele);
        point_ele = [point_ele(:),zeros(length_ele,1)];
    end
%% Transform boundary coords to boundary point coords.
    point_bnd = [bnd(1:2).', bnd(3)+zeros(2,1);
                 bnd(1:2).', bnd(4)+zeros(2,1)];
    point_bnd(:,3) = 0; % z coordinate
    point_bnd(:,4) = 5000; % meshfactor  distance/2
	n_point_bnd = size(point_bnd, 1);
%% Checker Points
    if isempty(start_x)
        if le(max(size(sizes)),2)
            start_x = -anz(end)*sizes(2)*0.5;
        else
            for it = 1: anz(2)
                start_x = sizes(anz(1)+it);
            end
            start_x = -start_x*0.5;
        end
    end

    y_check = zeros(1,anz(1));
    x_check = zeros(1,anz(end));
    end_x = start_x+sum(sizes(anz(1)+1:end));
    if le(max(size(sizes)),2)
        y_check = [0,cumsum(ones(1,anz(1))*sizes(1))];
        x_check = [start_x, start_x+cumsum(ones(1,anz(end))*sizes(end))];

    else
        if sum(anz)~=length(sizes)
            error('Sizes and numbers of checkers do not fit together.');
        end
        y_check = [0,cumsum(sizes(1:anz(1)))];
        x_check = [start_x, start_x+cumsum(sizes(anz(1)+1:end))];
    end


    point_checkers = [];
    mf_point_checkers = [];
    if ~le(max(size(sizes)),2)
        %temporary additional values to  make mesh factor computation easier
        sizes_temp = [10000,sizes(1:anz(1)),10000,10000,sizes((anz(1)+1):sum(anz)),10000];
        anz_temp=anz+2; %adjusted for sizes_temp
    end
    for it = 1: length(x_check) % mesh factor computation
        for itt = 1: length(y_check)
            point_checkers = [point_checkers;...
                x_check(it),y_check(itt)];
            if le(max(size(sizes)),2)
                mf_point_checkers = [mf_point_checkers;...
                    min(sizes(1),sizes(end))]; %#ok<*AGROW>
            else
                mf_point_checkers = [mf_point_checkers;...
                    min(min(sizes_temp(anz_temp(1)+it),sizes_temp(anz_temp(1)+it+1)),min(sizes_temp(itt),sizes_temp(itt+1)))];
               %
            end
        end
    end

%% Remove Ckecker points points which coincide with electrode positions.
    %ele_at_checker = all(ismember(point_checkers, point_ele).').'; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%no
    [ele_at_checker,~] = ismember(point_checkers,point_ele,'rows');
    point_checkers(ele_at_checker, :) = [];

    n_point_checkers = size(point_checkers, 1);
    n_point_ele = size(point_ele, 1);

    mf_point_ele = ones(1,n_point_ele)*0.05;
    mf_point_checkers(ele_at_checker, :)= [];

    point_ele(:,3:4) = [zeros(length_ele,1), ones(length_ele,1)*0.05]; % z Koordinate, mashfaktor
    point_checkers(:,3:4)= [zeros(n_point_checkers,1), mf_point_checkers];


 %% Create point list
    point_bnd   = [(1:n_point_bnd)', point_bnd];
    point_ele   = [(n_point_bnd+1:(n_point_bnd+n_point_ele))', point_ele];
    point_checkers = [((n_point_bnd+n_point_ele+1):(n_point_bnd+n_point_ele+n_point_checkers))', point_checkers];
    points = [point_bnd; point_ele; point_checkers];

    a = eq(points(:,3),0);
    b = ge(points(:,2),-60);
    c = le(points(:,2),60);
    mesh_correction = logical(a.*b.*c);
    %if checker points lie between receivers the mesh factor is lowered as
    %to not mess with the mesh generator (it gets upset if points with
    %widely different mesh factors lie to close together
    points(mesh_correction,5) = 0.05;
    clear a b c
%% Create lines
    %all lines point from left to right or downwards
    lines = [];

    %B = sortrows(A,column) sorts rows in ascending order based on column
    lines_bnd = [1 ,3; 3,4;2,4]; %left, down, right


    temp_points =  points(eq(points(:,3),0),:);
    [~,x_sort] = sortrows(temp_points,2);

    %lines upper boundary
    lines_bnd = [lines_bnd; ...
                temp_points(x_sort(1:end-1),1),temp_points(x_sort(2:end),1)];
    lines_bnd = [(1:size(lines_bnd,1))', lines_bnd];
    lines_count = size(lines_bnd,1)+1;

    clear x_sort temp_points;
    checker_pts_no_bnd = points(5:end,:); %remove corner points
    checker_pts_no_bnd =checker_pts_no_bnd(find(checker_pts_no_bnd(:,3)>0),:); %#ok<FNDSB>
    %all checker pts not on boundary

    %horizontal checker lines
    lines_checker_hor=[];
    [~,y_sort] = sortrows(checker_pts_no_bnd,3);
    for i = 1:(length(y_sort)/(anz(2)+1)) %Lines per depth
        st = 1+(i-1)*(anz(2)+1);
        lines_checker_hor= [lines_checker_hor; ...
                        checker_pts_no_bnd(y_sort(st:(st+anz(2)-1)),1), ...
                        checker_pts_no_bnd(y_sort((st+1):(st+anz(2))),1) ...
            ];
    end
    lines_checker_hor=[(lines_count:(lines_count+size(lines_checker_hor,1)-1))',lines_checker_hor];
    lines_count = lines_count+size(lines_checker_hor,1);
    %vertical checker lines


    [tmp_checker_in_ele,~] = ismember(points(:,2),x_check);
    checker_pts_bnd = points(tmp_checker_in_ele,:);
    checker_pts_bnd = checker_pts_bnd(eq(checker_pts_bnd(:,3),0),:); %checker points at boundary
    points_checker_raw = [checker_pts_bnd; checker_pts_no_bnd];
    [~,x_sort] = sortrows(points_checker_raw,2);
    lines_checker_vert = [];
    for i = 1:(length(x_sort)/(anz(1)+1)) %Lines per x
        st = 1+(i-1)*(anz(1)+1);
        lines_checker_vert= [lines_checker_vert; ...
                        points_checker_raw(x_sort(st:(st+anz(1)-1)),1), ...
                        points_checker_raw(x_sort((st+1):(st+anz(1))),1) ...
            ];
    end
    lines_checker_vert = [(lines_count:(lines_count+size(lines_checker_vert,1)-1))',lines_checker_vert];

    lines = [lines_bnd;lines_checker_hor; lines_checker_vert];

    %lines = [(1:size(lines,1))', lines]; %schon fr�her index
    %lines_bnd = [(1:size(lines_bnd,1))', lines_bnd];
%% lower surface
    surfaces = [];
    %surfaces = sparse(surfaces);
%     first surface
%     _______     ________x (index =2)
%     |      |___|        |
%     |                   |
%     |___________________|
    big_surf = [1,2,-3];
    point_A = 2;            %startng point

    % all lines to the right of the checker field
    while ~eq(x_check(end),points(point_A,2))
        [point_B_row,~] = ismember(lines_bnd(:,3), point_A);
        if ~isscalar(lines_bnd(point_B_row,2))
                error('Too many matching points.');
        else
            big_surf = [big_surf, -lines_bnd(point_B_row,1)];
            point_A = lines_bnd(point_B_row,2);
        end
    end
    %right dive
    % all vetical lines of the right of the checker field
    while ~eq(y_check(end),points(point_A,3))
        [point_B_row,~] = ismember(lines_checker_vert(:,2), point_A);
        if ~isscalar(lines_checker_vert(point_B_row,2))
            error('Too many matching points.');
        else
            big_surf = [big_surf, lines_checker_vert(point_B_row,1)];
            point_A = lines_checker_vert(point_B_row,3);
        end
    end
    %lower checker bnd
    while ~eq(x_check(1),points(point_A,2))
        [point_B_row,~] = ismember(lines_checker_hor(:,3), point_A);
        if ~isscalar(lines_checker_hor(point_B_row,3))
                error('Too many matching points.');
        else
            big_surf = [big_surf, -lines_checker_hor(point_B_row,1)];
            point_A = lines_checker_hor(point_B_row,2);
        end
    end
    %left checker
    while ~eq(0,points(point_A,3))
        [point_B_row,~] = ismember(lines_checker_vert(:,3), point_A);
        if ~isscalar(lines_checker_vert(point_B_row,3))
            error('Too many matching points.');
        else
            big_surf = [big_surf, -lines_checker_vert(point_B_row,1)];
            point_A = lines_checker_vert(point_B_row,2);
        end
    end
    %left of checker
    while ~eq(1,point_A)
        [point_B_row,~] = ismember(lines_bnd(:,3), point_A);
        if ~isscalar(lines_bnd(point_B_row,2))
                error('Too many matching points.');
        else
            big_surf = [big_surf, -lines_bnd(point_B_row,1)];
            point_A = lines_bnd(point_B_row,2);
        end
    end

%% checker surfaces

    checker_surf = zeros((anz(1)*anz(2)),(length_ele+4)); % 121+4 = maximal amount of lines possible in a checker
    surf_row = 1;


    for yy = y_check(1:(end-1))
        for xx = x_check(1:(end-1))
            surf_col = 1;
            [point_temp1_row,~] = ismember(points(:,2), xx);
            [point_temp2_row,~] = ismember(points(:,3), yy);
             point_A_row = logical(point_temp1_row.*point_temp2_row); %determine upper left point of current checker field
             point_A = points(point_A_row,1);
             point_A_control = point_A;
            %left down
            [point_B_row,~] = ismember(lines_checker_vert(:,2), point_A);
            if ~isscalar(lines_checker_vert(point_B_row,3))
                error('Too many matching points 1.');
            else
                checker_surf(surf_row,surf_col) = lines_checker_vert(point_B_row,1);
                surf_col = surf_col+1;
                point_A = lines_checker_vert(point_B_row,3);
            end

            %low ->
            [point_B_row,~] = ismember(lines_checker_hor(:,2), point_A);
            if ~isscalar(lines_checker_hor(point_B_row,3))
                    error('Too many matching points 2.');
            else
                checker_surf(surf_row,surf_col) = lines_checker_hor(point_B_row,1);
                surf_col = surf_col+1;
                point_A = lines_checker_hor(point_B_row,3);
            end

            %right up
            [point_B_row,~] = ismember(lines_checker_vert(:,3), point_A);
            if ~isscalar(lines_checker_vert(point_B_row,2))
                error('Too many matching points 3.');
            else
                checker_surf(surf_row,surf_col) = -lines_checker_vert(point_B_row,1);
                surf_col = surf_col+1;
                point_A = lines_checker_vert(point_B_row,2);
            end

            %top <-
            if eq(yy,0)
                while ~eq(point_A_control,point_A)
                    [point_B_row,~] = ismember(lines_bnd(:,3), point_A);
                    if ~isscalar(lines_bnd(point_B_row,2))
                            error('Too many matching points 4.');
                    else
                        checker_surf(surf_row,surf_col) = -lines_bnd(point_B_row,1);
                        surf_col = surf_col+1;
                        point_A = lines_bnd(point_B_row,2);
                    end
                end
            else
                [point_B_row,~] = ismember(lines_checker_hor(:,3), point_A);
                if ~isscalar(lines_checker_hor(point_B_row,2))
                        error('Too many matching points 4.');
                else
                    checker_surf(surf_row,surf_col) = -lines_checker_hor(point_B_row,1);
                    surf_col = surf_col+1;
                    point_A = lines_checker_hor(point_B_row,2);
                end

            end
            surf_row = surf_row+1;
        end
    end

%% Define Physical lines
    point_A = 1;
    y_min = [];
    lines_bnd_top = lines_bnd(4:end,:);
    while ~eq(point_A,2)
        [point_B_row,~] = ismember(lines_bnd_top(:,2), point_A);
        if ~isscalar(lines_bnd_top(point_B_row,3))
                error('Too many matching points 5.');
        else
            y_min = [y_min, lines_bnd_top(point_B_row,1)];
            point_A = lines_bnd_top(point_B_row,3);
        end
    end



    y_max = big_surf(2);
    x_min = big_surf(1);
    x_max = -big_surf(3);


%% write Geo file
    fileID = fopen([gmsh_file, '.geo'], 'w');

    % Handle geometric entities.
    % Note: multiple defined lines are merged automatically in Gmsh!

    % Add points.
    for ii = 1:size(points,1)
        fprintf(fileID, 'Point(%d) = {%d, %d, %d, %d};\n', ...
            points(ii,:));
    end
    % Add lines
    fprintf(fileID, '\n');
    for ii = 1:size(lines, 1)
        fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
            lines(ii,:));
    end

    % Add Line loops
    fprintf(fileID, 'Line Loop(1) = ');
    fprintf(fileID,strrep(['{' sprintf(' %d,', big_surf) ')'], ',)', '};'));
    fprintf(fileID, '\n');

    for ii = 1:size(checker_surf,1)
        temp_surf = nonzeros(checker_surf(ii,:));
        fprintf(fileID, 'Line Loop(%d) = ', (ii+1));
        fprintf(fileID,strrep(['{' sprintf(' %d,', temp_surf) ')'], ',)', '};'));
        fprintf(fileID, '\n');
    end
    for i = 1:(1+size(checker_surf,1))
        fprintf(fileID,'Plane Surface(%d) = {%d};\n', i, i);
    end

    fprintf(fileID, 'Physical Line("xmin") = ');
    fprintf(fileID,strrep(['{' sprintf(' %d,', x_min) ')'], ',)', '};')) ;
    fprintf(fileID, '\n');
    fprintf(fileID, 'Physical Line("xmax") = ');
    fprintf(fileID,strrep(['{' sprintf(' %d,', x_max) ')'], ',)', '};')) ;
    fprintf(fileID, '\n');
    fprintf(fileID, 'Physical Line("ymin") = ');
    fprintf(fileID,strrep(['{' sprintf(' %d,', y_min) ')'], ',)', '};')) ;
    fprintf(fileID, '\n');
    fprintf(fileID, 'Physical Line("ymax") = ');
    fprintf(fileID,strrep(['{' sprintf(' %d,', y_max) ')'], ',)', '};')) ;
    fprintf(fileID, '\n');

   dom_name1 = 'res_1';
   dom_name2 = 'res_2';

   if all(group == 0)
        fprintf(fileID, ...
        'Physical Surface("%s") = ',dom_name1);
        fprintf(fileID,strrep(['{' sprintf(' %d,', 1:(1+size(checker_surf,1))) ')'], ',)', '};'));
   elseif all(group == 1)
        fprintf(fileID, ...
        'Physical Surface("%s") = ',dom_name1);
        fprintf(fileID,'{1};\n');
        % big_surf is always res_1
        fprintf(fileID, ...
        'Physical Surface("%s") = ',dom_name2);
        fprintf(fileID,strrep(['{' sprintf(' %d,', 2:(1+size(checker_surf,1))) ')'], ',)', '};'));
   else
       group = reshape(group.',1,[]);
       res_1_index = find(group==0);
       res_2_index = find(group==1);

       res_1_index = [1, (res_1_index +1)]; % um index big_surf zu ber�cksichtigen
       res_2_index = res_2_index +1;

       fprintf(fileID, ...
        'Physical Surface("%s") = ',dom_name1);
       fprintf(fileID,strrep(['{' sprintf(' %d,',res_1_index) ')'], ',)', '};')) ;
       fprintf(fileID, '\n');
       fprintf(fileID, ...
        'Physical Surface("%s") = ',dom_name2);
       fprintf(fileID,strrep(['{' sprintf(' %d,', res_2_index) ')'], ',)', '};'));
   end

   fclose(fileID);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Gmsh.

    system([gmsh_path.folder, '/gmsh -2 ', ...
            gmsh_file, '.geo -v 0 -format msh2']);

% gmsh -2 sens5.geo -v 0 -format msh2
% gmsh -2 .geo -v 0 -format msh2
    %% Import mesh information from .msh file.

    mesh = Mesh.loadGmsh([gmsh_file, '.msh']);

    % Override type.
    mesh.type = 'gmsh_load';

    %% Clean up.

%    delete([gmsh_file, '.geo'], [gmsh_file, '.msh']);

end


