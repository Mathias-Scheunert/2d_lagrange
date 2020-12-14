function A = getCellArea(mesh)
    % Calculates the triangle areas of cells in mesh using Heron's formula.
    %
    % SYNTAX
    %    a = getCellArea(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OUTPUT PARAMETER
    %   a ... Vector, of mesh cell areas.

    assert(isfield(mesh, 'cell2cord'));
    % TODO: If not available, just calculate.

    % Fetch.
    n_cell = size(mesh.cell2vtx, 1);

    % Loop over cells and calc areas.
    A = zeros(n_cell, 1);
    for ii = 1:n_cell

        % Fetch points.
        p1 = mesh.cell2cord{ii}(1, :);
        p2 = mesh.cell2cord{ii}(2, :);
        p3 = mesh.cell2cord{ii}(3, :);

        % Fetch edges.
        a = norm(p1 - p2);
        b = norm(p1 - p3);
        c = norm(p2 - p3);
        s = (a+b+c)/2;

        % Get area.
        A(ii) = sqrt(s*(s - a)*(s - b)*(s - c));
    end
end
