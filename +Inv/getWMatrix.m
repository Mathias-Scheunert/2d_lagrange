function W = getWMatrix(mesh)
    % Calculates the representation of the smoothness operator in FE space.
    %
    % Derivative follows C. Schwarzbach, Finite element based inversion for
    %   time-harmonic electromagnetic problems, 2013
    %
    % SYNTAX
    %   W = getWMatrix(mesh)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OUTPUT PARAMETER
    %   W ... Matrix representation of the smoothness operator acting on 
    %         the vector of model parameters, i.e. the weak form of the 
    %         Poisson problem in the space of the piece-wise constant
    %         parameter (function), represented by a mixed formulation
    %         (Raviart-Thomas elements).
   
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx'})), ...
        'mesh - appended struct, containing cells info, expected.');
    
    %% Set up matricies for Raviart-Thomas element Poisson problem.
    
    % Initialize RT finite elemetns.
    feRT = FeRT.initFiniteElement(mesh);
    
    % Assemble mass and divergence operator matrices.
    M = FeRT.assembleMass(feRT, mesh);
    D = FeRT.assembleDiv(feRT, mesh);
    
    %% Derive weighting matrix.
    
    R = chol(M);
    W = R \ D.';
end

