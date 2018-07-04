function bnd = assignBC(bnd, fe, mesh, verbosity)
    % Assigns (or loads) the values to each boundary DOF.
    %
    % SYNTAX
    %   bnd = assignBC(bnd, fe, mesh, verbosity)
    %
    % INPUT PARAMETER
    %   bnd  ... Struct, containing the preset of the BC types.
    %   fe   ... Struct, including all information to set up Lagrange FE,
    %            as well as the linear system components.
    %   mesh ... Struct, containing mesh information, i.e. coordinates
    %            of vertices and its relation to the triangles and edges.
    %
    % OPTIONAL PARAMETER
    %   verbosity ... Logical, denoting if current status should be
    %                 printed.
    %
    % OUTPUT PARAMETER
    %   bnd ... Struct, containing complete boundary DOF information.
    %
    % REMARKS
    %   Three types of given bnd information for the sides of a predefined
    %   mesh can be handled:
    %       scalar          - consant value is set on every edge DOF
    %       vector          - each element is assigned to the respective 
    %                         edge DOF
    %       function handle - the value is calculated from the @fun... for 
    %                         (the coordinate of) each edge DOF

    %% Check input.
    
    assert(isstruct(fe) && all(isfield(fe, {'order', 'sizes'})), ...
        'fe - struct, including all information of FE linear system, expected.');
    assert(isstruct(mesh) && all(isfield(mesh, {'vertices', 'edge2vtx', 'bnd_edge_bot'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isstruct(bnd) && all(isfield(bnd, {'type', 'val'})), ...
        'bnd - struct, containing basic boundary condition information, expected.');
    if nargin < 4
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end
    
    %% Get or load bnd DOF.
    
    if verbosity
       fprintf('Get or load bnd DOFs ... '); 
    end
    
    % Get bndDOF.
    switch mesh.type
        case {'cube', 'rhomb'}          
            if ~isfield(bnd, 'bndDOF')
                bnd.bndDOF = Fe.getBndDOF(fe, mesh);
            end
            
        case {'external'}
            error('Not implemented yet.');
            % TODO: if mesh is given, load information at this stage.
            
        otherwise
            error('BC handling only supported for kown mesh type.');
    end
    if verbosity
       fprintf('done.\n'); 
    end
    
    %% Assign values for bnd DOFs.
    
    if verbosity
       fprintf('Assign or load bnd DOF values ... '); 
    end
    
    switch mesh.type
        case {'cube', 'rhomb'}
            
            % Check if every domain boundary will be handled.
            check_BC = cell2mat(cellfun(@(cur_val) ...
                {cellfun(@(x) ~isempty(x), cur_val)}, bnd.val));
            check_BC = sum(check_BC, 2);
            assert(all(check_BC), ...
                'Missing BC for domain boundaries. Please check DRIVE.');
            
            % Iterate over the different bnd types.
            for ii = 1:length(bnd.type)
                % Notice empty info.
                empty_BC = cellfun(@isempty, bnd.val{ii});
                
                % Notice info, given in basic (i.e. scalar form).
                basic_BC = cellfun(@(x) ...
                    isscalar(x) && ~isa(x, 'function_handle'), ...
                    bnd.val{ii});
                
                % Notice info, given as function handle.
                fun_BC = cellfun(@(x) isa(x, 'function_handle'), ...
                    bnd.val{ii});

                % Calculate Dirichlet values from given function handle.
                if any(fun_BC)
                    bnd.val{ii}(fun_BC) = cellfun(@(f, coo) ...
                        {arrayfun(@(x, y) f(x, y), coo(:, 1), coo(:, 2))}, ...
                        bnd.val{ii}, bnd.bndDOF.bnd_DOF_coo);
                end
                                
                % Replace basic info with bndDOF related values.
                if any(basic_BC)
                    bnd.val{ii}(basic_BC) = cellfun(@(x, y) {x + zeros(size(y))}, ...
                        {bnd.val{ii}{basic_BC}}, {bnd.bndDOF.bnd_DOF{basic_BC}});
                end
                
                % Check, if given bndDOF related values have correct size.
                DOF_related_BC = ~empty_BC & ~basic_BC & ~fun_BC;
                if any(DOF_related_BC)
                    consistent_BC = cellfun(@(x, y) length(x) == length(y), ...
                        {bnd.val{ii}{DOF_related_BC}}, ...
                        {bnd.bndDOF.bnd_DOF{DOF_related_BC}}).';
                    assert(isempty(consistent_BC) || all(consistent_BC), ...
                        ['Number of given DOF related BC values do not ', ...
                        'match the bnd_DOF.']);
                end
            end
            
        case {'external'}
            % TODO: implement.
    end
    if verbosity
       fprintf('done.\n'); 
    end
end