function fe = initFiniteElement(order, mesh, point, verbosity)
    % Provide finite element structure for 2D Lagrange elements.
    %
    % SYNTAX
    %   fe = initFiniteElement(order, mesh, point[, verbosity])
    %
    % INPUT PARAMETER
    %   order ... Scalar, denoting the order of elements.
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %   point ... Matrix nx2, containing coordinates of observation points.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   fe ... Struct, including all information to set up Lagrange FE.
    %
    % REMARKS
    %   Please note the ordering of vertex and edge indices within
    %   Fe.getDOFMap.m!
    
    %% Check input
    
    assert(isscalar(order), ...
        'Imput of scalar order expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx'})), ...
        'mesh - appended struct, containing cells info, expected.');
    if mesh.dim ~= 2
        error('Only 2D problem supported yet.');
    end
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    %% Collect information.
    
    fe = struct();
    fe.dim = mesh.dim;
    fe.order = order;

    % Get quadrature rules.
    if verbosity
       fprintf('Set up quadrature rule ... '); 
    end
    [fe.quad.nodes, fe.quad.weights] = Quad.getQuadratureRule(fe.order, fe.dim);
    if verbosity
       fprintf('done.\n'); 
    end
    
    % Get Lagrange basis functions and its gradients.
    if verbosity
       fprintf('Set up Lagrange basis ... '); 
    end
    fe.base = FeL.getBasis(fe.order);
    if verbosity
       fprintf('done.\n'); 
    end
    
    % Get DOF mapping.
    if verbosity
       fprintf('Obtain DOF maps ... '); 
    end
    fe.DOF_maps = FeL.getDOFMap(mesh, fe);
    if verbosity
       fprintf('done.\n'); 
    end
    
    % Get common sizes.
    fe.sizes.quad_point = size(fe.quad.nodes, 1);
    fe.sizes.cell = size(mesh.cell2vtx, 1);
    fe.sizes.vtx = size(mesh.vertices, 1);
    fe.sizes.DOF = fe.DOF_maps.n_DOF;
    fe.sizes.DOF_loc = size(fe.base.DOF, 1);
    
    % Get interpolation operator.
    if verbosity
       fprintf('Obtain interpolation operator ... '); 
    end
    fe.I = FeL.getInterpolation(fe, mesh, point);
    if verbosity
       fprintf('done.\n');
       fprintf('... FE struct initialized.\n \n');
    end
end