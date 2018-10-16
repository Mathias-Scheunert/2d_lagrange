function D = assembleDiv(fe, mesh, verbosity)
    % Assembles the sparse divergence operator matrix.
    %
    % continuous:
    %  f(v) = \int_Omega \div u d(x,y) = ...
    % Galerkin approx. 
    % (i.e. piece-wise evaluation w.r.t. simplices including the  
    %     numerical quadrature approx for integral evaluation and coord.
    %     shift to reference simplex):
    %   \sum_k \sum_l ( w_l q_k div(\phi_i({x,y}_l)) 1/\abs(\det(B_k)) ) 
    %   * \abs(\det(B_k))
    % with
    %   q_k = 1 for all parameters in the triangle k and
    %   q_k = 0 elsewhere
    %
    % k - num simplices
    % l - num quadrature nodes
    % i - num basis functions    
    %
    % Using elemente-wise procedure to set up the global divergence matrix.
    %
    % SYNTAX
    %   D = assembleDiv(fe, mesh, param, verbosity)
    %
    % INPUT PARAMETER
    %   fe   ... Struct, including all information to set up Lagrange FE.
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   D  ... Matrix, representing the divergence part of the
    %          variational formulation.
    %
    % REMARKS
    %   As the basis (q_k) of the parameter space is 1 for all points
    %   within a triangle, i.e. for all quadrature nodes used within 
    %   assembling loop, it is not incorporated in the code below.
    % TODO: Implement in order to have a accurate/consistent description.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base'})), ...
        'fe - struct, including Lagrange reference element info , expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord', 'maps'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    if nargin < 3
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end

    %% Assemble divergence matrix.

    if verbosity
       fprintf('Assemble divergence matrix ... '); 
    end
       
    % Get common sizes.
    n_cell = fe.sizes.cell;
    n_DOF_glob = fe.sizes.DOF;
    n_entry_loc = fe.sizes.DOF_loc;
    
    % Extract quadrature info from struct.
    gauss_cords = fe.quad.nodes;
    gauss_weights = num2cell(fe.quad.weights.');
        
    % Initialize index and value vector for sparse matrix assembling.
    [j, i, s] = deal(zeros(n_cell * n_entry_loc, 1));
    
    % Set up recurring quantity.
    % Get divergence of basis functions for all quadrature nodes referred 
    % to the reference simplex.
    basis_eval = arrayfun(@(x,y) {fe.base.div_Phi(x, y)}, ...
                    gauss_cords(:,1), gauss_cords(:,2)).';
                
    % Iterate over all simplices.
    for ii = 1:n_cell
        % W.r.t. global coords:
        % To achive normal component of basis functions of adjacent cells
        % to be consistent, check if the related edge normal is parallel 
        % (or antiparallel) oriented with the global normal of the edge
        % (see FeRt.initFiniteElement for its definition).
        %     parallel: Nothing to do. 
        % antiparallel: Switch sign of the basis function.
        % Note: At first, normals are ordered w.r.t. to the mesh.cell2edge 
        % ordering (referred to the relations of cells and edges in global 
        % coords)
        % Get all edges global normal of current cell.
        cur_edge_global_normals = fe.glo_edge_normals(mesh.cell2edg(ii,:), :);
        % Get current edge local normal.
        cur_edge_normals = Mesh.getEdgeNormal(mesh, ...
                               mesh.cell2edg(ii,:), ii);
        cur_edge_normals = cell2mat([cur_edge_normals{:}].');
        % Obtain sign function.
        Phi_sign = dot(cur_edge_global_normals, cur_edge_normals, 2);
        % Change ordering such that sign function can be applied on (local)
        % basis functions for reference simplex.
        Phi_sign = Phi_sign(mesh.loc2glo).';
        
        % Set up kernel for integral (quadrature summation).
        % By multiplying the vector of basis functions with the
        % reference/source function.
%         quad_kern = cellfun(@(x, y) {x * y * 1/(mesh.maps{ii}.detB)}, ... 
%                         basis_eval, gauss_weights); 
        quad_kern = cellfun(@(x, y) ...
                        {(Phi_sign.*x) * y * 1/abs(mesh.maps{ii}.detB)}, ...
                        basis_eval, gauss_weights);
        
        % Evaluate numerical integration and incorporate the norm of the 
        % Jacobi-determinat due to mapping back from reference simplex to 
        % global coordinates.
        m_loc = abs(mesh.maps{ii}.detB) * sum(cat(3, quad_kern{:}), 3);
                  
        % Fill up index and value vectors.
        j_loc = fe.DOF_maps.cell2DOF{ii};
        glob_idx_start = ((ii - 1) * n_entry_loc) + 1;
        glob_idx_end = glob_idx_start + n_entry_loc - 1;
        j(glob_idx_start:glob_idx_end) = j_loc(:);
        i(glob_idx_start:glob_idx_end) = ii + (0 * j_loc(:));
        s(glob_idx_start:glob_idx_end) = m_loc(:);
    end
    
    % Create sparse rhs vector.
    D = sparse(i, j, s, n_cell, n_DOF_glob);
    
    if verbosity
       fprintf('done.\n'); 
    end
end