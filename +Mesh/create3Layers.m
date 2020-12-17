function mesh = create3Layers(ele,layers_y,verbosity, ref)
% Stores mesh information of a mesh created by Gmsh.
    %
    % Supported (tested) Gmsh version: 3.x, MSH file format version 2
    %
    % Therefore,
    %   1) a .geo input file is created
    %   2) a .msh output file is obtained from running Gmsh
    %   3) information from .msh is converted
    %   4) .geo and .msh files are deleted
    %
    % SYNTAX
    %   mesh = createGmsh(bnd, args)
    %
    % INPUT PARAMETER
    %   ele  ... electrode positions nx2
    %   layers_y. starting depths of 2nd and 3rd layer must be positive and
    %             smaller than bnd values
    % OUTPUT PARAMETER
    %   mesh ... Struct, containing vertex2coordinates (vertices),
    %            simplex2vertex (cell2vtx), boundary egde2coordinates
    %            (bnd_edge2vtx), boundary egde 2 model domain boundary
    %            (bnd_edge_...) and simplex parameter domains
    %            (parameter_domains).
    %            domain_names =  res_1, res_2, res_3

    %% Check Gmsh version.

    gmsh_path = dir('**/gmsh.exe');
    if isempty(gmsh_path)
        error('Gmsh executable could not be found in path.');
    end

    %% Create Gmsh input file.

    gmsh_file = 'three_layers';
    createGmshInput([gmsh_file, '.geo'], ele, layers_y);

    %% Run Gmsh.

     system([gmsh_path.folder, '/gmsh -2 ', ...
            gmsh_file, '.geo -v 0 -format msh2 ']);

    %% Import mesh information from .msh file.

    mesh_type = 'gmsh_load';
    mesh = Mesh.initMesh(mesh_type, 'name', [gmsh_file, '.msh'], ...
                 'ref', ref, 'verbosity', verbosity);

    delete([gmsh_file, '.geo'], [gmsh_file, '.msh']);
end

function createGmshInput(name, ele, layers_y)
 % Creates an .geo input file for Gmesh
 %
 % REMARKS
 %  Note that the 2D grid will be referred to gmesh x and z coordinate:
 %  0,0 ------>
 %      |      x = gmesh-x
 %      |
 %      |
 %      v
 %       y = gmesh-z
 %
%% Set domain boundaries.

TXRX = unique(ele, 'rows');
if ~isempty(TXRX)

    % Get specific extentions.
%    min_x = min(TXRX(:,1));
%    max_x = max(TXRX(:,1));
    offset = max(TXRX(:,1))-min(TXRX(:,1));

    % Adapt.

    bnd = [-600*offset, ...
           600*offset, ...
           0, ...
           600*offset];
end

% Transform boundary coordinates.
X = [bnd(1:2), fliplr(bnd(1:2))];
Y = [bnd(3:4); bnd(3:4)];
Y = Y(:).';

%% Set up domain basic geometry entities (points, lines).
% Create point input ([n x 4] matrix) structure for domain boundaries.
n_domain = length(X);
points = zeros(n_domain, 4);
for i = 1:n_domain
   points(i, 1) = i ;       % point number
   points(i, 2) = X(i);     % gmesh-x cooridnate.
   points(i, 3) = 0;        % gmesh-y cooridnate.
   points(i, 4) = Y(i);     % gmesh-z cooridnate.
end
% Add Layer depths
points(n_domain+1, :) = [n_domain+1,bnd(1),0,layers_y(1)];
points(n_domain+2, :) = [n_domain+2,bnd(1),0,layers_y(2)];
points(n_domain+3, :) = [n_domain+3,bnd(2),0,layers_y(1)];
points(n_domain+4, :) = [n_domain+4,bnd(2),0,layers_y(2)];

n_domain = n_domain+4;
% Add TX/RX positions.
for i = 1:size(ele,1)
   points(i+n_domain, 1) = i+n_domain ;   % point number
   points(i+n_domain, 2) = ele(i,1);      % gmesh-x cooridnate.
   points(i+n_domain, 3) = 0;             % gmesh-y cooridnate.
   points(i+n_domain, 4) = ele(i,2);      % gmesh-z cooridnate.
end
% Create line input ([n x 3] matrix)
points_temp = sortrows(points,[4 2 1]);
temp_idx = points_temp((points_temp(:,4)==0),1); %point idx on surface

lines_hor = [temp_idx(1:end-1),temp_idx(2:end)];
clear temp_idx
temp_idx = points_temp((points_temp(:,4)==layers_y(1)),1); %point idx second layer
lines_hor = [lines_hor; ...
            temp_idx(1:end-1),temp_idx(2:end)];
clear temp_idx
temp_idx = points_temp((points_temp(:,4)==layers_y(2)),1); %point idx third layer
lines_hor = [lines_hor; ...
            temp_idx(1:end-1),temp_idx(2:end)];
clear temp_idx
temp_idx = points_temp((points_temp(:,4)==bnd(4)),1); %point idx bottom
lines_hor = [lines_hor; ...
            temp_idx(1:end-1),temp_idx(2:end)];
clear temp_idx  points_temp

points_temp = sortrows(points,[2 4 1]);
temp_idx = points_temp((points_temp(:,2)==bnd(1)),1); %point left
lines_vert = [temp_idx(1:end-1),temp_idx(2:end)];
clear temp_idx

temp_idx = points_temp((points_temp(:,2)==bnd(2)),1); %point left
lines_vert = [lines_vert; ...
              temp_idx(1:end-1),temp_idx(2:end)];
clear temp_idx points_temp
lines = [lines_hor;lines_vert];
lines = [(1:size(lines,1))',lines];

% Create Surfaces
surf = cell(3,1);
surf{1} = [1:(size(ele,1)+1), lines((size(lines_hor,1)+4),1),-lines((size(ele,1)+2),1),-lines((size(lines_hor,1)+1),1)];
surf{2} = [lines((size(ele,1)+2),1),lines((size(lines_hor,1)+5),1),-lines((size(ele,1)+3),1),-lines((size(lines_hor,1)+2),1)];
surf{3} = [lines((size(ele,1)+3),1),lines((size(lines_hor,1)+6),1),-lines((size(ele,1)+4),1),-lines((size(lines_hor,1)+3),1)];

%% Create Gmsh geometry (text) file and define higher order entities.

fileID = fopen(name, 'w');
fprintf(fileID, '// definitions \n');
fprintf(fileID, 'max_off = %d;\n',offset);
dom_wdt = max(bnd(2)-bnd(1));
fprintf(fileID, 'dom_wdt = %d;\n',dom_wdt);
%for distance fields
fprintf(fileID, 'lc_min = 0.05;\nlc_max = dom_wdt/2;\nDistMax = dom_wdt/0.5;\nDistMin = 0.8;\n Mesh.Algorithm=1; \n');

% Add geometry-describing point ids.
for i = 1:2
    fprintf(fileID, 'Point(%d) = {%d, %d, %d};\n', ...
        points(i,:));
end
for i = 3:4
    fprintf(fileID, 'Point(%d) = {%d, %d, %d, dom_wdt/2};\n', ...
        points(i,:));
end
for i = 5:size(points,1)
    fprintf(fileID, 'Point(%d) = {%d, %d, %d};\n', ...
        points(i,:));
end
% Add line ids.
fprintf(fileID, '\n');
for i = 1:size(lines,1)
    fprintf(fileID, 'Line(%d) = {%d, %d};\n', ...
        lines(i,:));
end

% Define polygonal chain describing the complete domain boundary.
fprintf(fileID, '\n');
for i = 1:3
    fprintf(fileID, ...
    ['Line Loop(%d) = {', repmat('%d, ', 1, length(surf{i}(1:end-1))), '%d};\n'], ...
    i, surf{i}(1:end));
end
fprintf(fileID, '\n');
fprintf(fileID, ...
        'Plane Surface(1) = {1};\nPlane Surface(2) = {2};\nPlane Surface(3) = {3};\n');
% Add boundary (physical) line ids.
% xmin (left) size(lines_hor,1)
fprintf(fileID, '\n');
fprintf(fileID,[ 'Physical Line("xmin") = {', repmat('%d, ', 1, 2), '%d};\n'], ...
     (size(lines_hor,1)+1):(size(lines_hor,1)+3));
fprintf(fileID, '\n');
fprintf(fileID,[ 'Physical Line("xmax") = {', repmat('%d, ', 1, 2), '%d};\n'], ...
     (size(lines_hor,1)+4):(size(lines_hor,1)+6));
fprintf(fileID, '\n');
fprintf(fileID, ...
    ['Physical Line("ymin") = {', repmat('%d, ', 1, (size(ele,1))), '%d};\n'], ...
     1:(size(ele,1)+1));
fprintf(fileID, '\n');
fprintf(fileID, 'Physical Line("ymax") = {%d};\n', size(ele,1)+4);


fprintf(fileID, '\n');
    fprintf(fileID, ...
        'Physical Surface("res_1") = {1};\nPhysical Surface("res_2") = {2};\nPhysical Surface("res_3") = {3};\n ');

% Add Distance Fields

fprintf(fileID, '\n');

ele_idx = n_domain+1 : (n_domain+size(ele,1));
for ii = 1:size(ele,1)
    fprintf(fileID,'Field[%d] = Distance;\n', 2*ii-1);
    fprintf(fileID,'Field[%d].NodesList = {%d};\n', 2*ii-1,ele_idx(ii));
    fprintf(fileID,'Field[%d] = Threshold;\n', 2*ii);
    fprintf(fileID,'Field[%d].IField = %d;\n', 2*ii, 2*ii-1);
    fprintf(fileID,'Field[%d].LcMin = lc_min;\n', 2*ii);
    fprintf(fileID,'Field[%d].LcMax = lc_max;\n', 2*ii);
    fprintf(fileID,'Field[%d].DistMin = DistMin;\n', 2*ii);
    fprintf(fileID,'Field[%d].DistMax = DistMax;\n \n', 2*ii);
end

fprintf(fileID,'Field[%d] = Min;\n', 2*length(ele_idx)+1);
fprintf(fileID,'Field[%d].FieldsList = {',2*length(ele_idx)+1);
fprintf(fileID, [ ...
        repmat('%d, ', 1, length(ele_idx)-1), '%d};\n'], ...
        2*(1:(length(ele_idx))));
fprintf(fileID,'Background Field = {%d};\n', 2*length(ele_idx)+1);

% Remove multiples.
fprintf(fileID, '\n');
fprintf(fileID, 'Coherence;');

fclose(fileID);
end






