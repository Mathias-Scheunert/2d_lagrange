function fe = initFiniteElement(mesh, verbosity)
    % Provide finite element structure for 2D Raviart-Thomas elements.
    %
    % SYNTAX
    %   fe = initFiniteElement(mesh, point[, verbosity])
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   fe ... Struct, including all information to set up Lagrange FE.
    %
    % REMARKS
    %   Please note the ordering of edge indices within Fe.getDOFMap.m!
    %   Basis function definitions from local to global coordinates
    %   (ommitting the scaling factor) and Piola transformation
    %   (preserving orthogonality of vectors) are obtained from:
    %       Computational Bases for RTk and BDMk on Triangles;
    %       Ervin, V.J.; 2012

    %% Check input

    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'dim', ...
                                     'edge2vtx', 'vertices'})), ...
        'mesh - appended struct, containing cells info, expected.');
    if mesh.dim ~= 2
        error('Only 2D problem supported yet.');
    end
    if nargin < 2
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end

    %% Collect information.

    fe = struct();
    fe.dim = mesh.dim;
    fe.order = 1;

    % Get quadrature rules.
    % Note: For lowest-order RT-elements, assembling of mass matrix
    % requires order == 2.
    if verbosity
       fprintf('Set up quadrature rule ... ');
    end
    [fe.quad.nodes, fe.quad.weights] = Quad.getQuadratureRule(2*fe.order, fe.dim);
    if verbosity
       fprintf('done.\n');
    end

    % Get RT basis functions and its divergence.
    if verbosity
       fprintf('Set up Raviart-Thomas basis ... ');
    end
    fe.base = FeRT.getBasis();
    if verbosity
       fprintf('done.\n');
    end

    % Get DOF mapping.
    if verbosity
       fprintf('Obtain DOF maps ... ');
    end
    fe.DOF_maps = FeRT.getDOFMap(mesh);
    if verbosity
       fprintf('done.\n');
    end

    % Get common sizes.
    fe.sizes.quad_point = size(fe.quad.nodes, 1);
    fe.sizes.cell = size(mesh.cell2vtx, 1);
    fe.sizes.vtx = size(mesh.vertices, 1);
    fe.sizes.DOF = fe.DOF_maps.n_DOF;
    fe.sizes.DOF_loc = size(fe.base.DOF, 1);
    if verbosity
       fprintf('done.\n');
       fprintf('... FE struct initialized.\n \n');
    end
end
