function b = assembleRHS(fe, mesh, TX)
    % Assembles the rhs vector for different source types.
    %
    % SYNTAX
    %   b = assembleRHS(fe, mesh, TX)
    %
    % INPUT PARAMETER
    %   fe ... Struct, including all information to set up Lagrange FE.
    %   mesh  ... Struct, containing mesh information, i.e. coordinates
    %             of vertices and its relation to the triangles and edges.
    %   TX ... Struct, containing the source information. 
    %
    % OUTPUT PARAMETER
    %   b ... Vector, representing the given source w.r.t. the DOF.
    
    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'base', 'maps', 'DOF_maps'})), ...
        'fe - struct, including all information to set up Lagrange FE, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2cord'})), ...
        'mesh - appended struct, containing cell2cord info, expected.');
    assert(isstruct(TX) && all(isfield(TX, {'val', 'type'})), ...
        'TX - struct, including all source information, expected.');
    
    %% Assemble rhs vector.
    
    switch TX.type
        case 'point'
            
            % Get (first) cell which belongs to TX point by MATLAB bultin.
            cell_idx = tsearchn(mesh.vertices, mesh.cell2vtx, TX.coo);
            if isempty(cell_idx)
                error('tsearchn failed to identify cell that belongs to TX');
                % TODO: try fix that by using the 'own' search like already
                % done within Fe.getInterpolation().
            end
            
            % Get affine map for respective cell w.r.t TX.
            map = Mesh.getAffineMap(cell_idx, mesh, TX.coo);

            % Evaluate basis functions at TX.
            base = fe.base.Phi(map.xy_ref(1), map.xy_ref(2)).';
            
            % Get DOF indices.
            DOFs = fe.DOF_maps.cell2DOF{cell_idx};
                      
            % Create rhs vector by weight respective DOFs with TX strength.
            b = sparse(DOFs, 1, TX.val * base, fe.sizes.DOF, 1);
            
        case 'homogen'
            
            % Create rhs vector.
            b = TX.val + zeros(fe.sizes.DOF, 1);
            
        case 'secondary'
            error('Only point sources are supported yet.');
            
        otherwise
            error('Unknown source type.');
    end
end