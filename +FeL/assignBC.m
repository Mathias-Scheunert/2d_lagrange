function bnd = assignBC(bnd, fe, mesh, param, verbosity)
    % Assigns (or loads) the values to each boundary DOF.
    %
    % SYNTAX
    %   bnd = assignBC(bnd, fe, mesh[, verbosity])
    %
    % INPUT PARAMETER
    %   bnd   ... Struct, containing the preset of the BC types.
    %   fe    ... Struct, including all information to set up Lagrange FE,
    %             as well as the linear system components.
    %   mesh  ... Struct, containing the mesh information.
    %             For a detailed description of the content of the mesh
    %             struct please read header of Mesh.initMesh.
    %   param ... Vector of constant cell parameter values.
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
    assert(isstruct(mesh) && all(isfield(mesh, {'vertices', 'edge2vtx', 'bnd_edge'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    assert(isvector(param) && length(param) == size(mesh.cell2vtx, 1), ...
        'param - Vector of constant cell parameter values, expected.');
    assert(isstruct(bnd) && all(isfield(bnd, {'type', 'val', 'name'})), ...
        'bnd - struct, containing basic boundary condition information, expected.');
    if nargin < 5
        verbosity = false;
    else
        assert(islogical(verbosity), ...
            'verbosity - logical, denoting if status should be printed, expected');
    end

    %% Get or load bnd DOF.

    if verbosity
       fprintf('Get or load bnd DOFs ... ');
    end

    % Appand bndDOF information.
    switch mesh.type
        case {'cube', 'disc', 'gmsh_create', 'gmsh_load'}
            if ~isfield(bnd, 'bndDOF')
                bnd = FeL.getBndDOF(fe, mesh, bnd);
            end
        otherwise
            error('BC handling only supported for known mesh types.');
    end
    if verbosity
       fprintf('done.\n');
    end

    %% Assign values for bnd DOFs.

    if verbosity
       fprintf('Assign or load bnd DOF values ... ');
    end

    % Check if every domain boundary will be handled.
    check_BC = cell2mat(cellfun(@(cur_val) ...
        {cellfun(@(x) ~isempty(x), cur_val)}, bnd.val));
    check_BC = sum(check_BC, 2);
    assert(all(check_BC), ...
        'Missing BC for domain boundaries. Please check DRIVE.');
    assert(all(check_BC < 2), ...
        'Some boundaries are assigned multiply - check DRIVE.');

    % Iterate over the different bnd types.
    for ii = 1:length(bnd.type)
        switch bnd.type{ii}
            case 'dirichlet'
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
            % Note: a scalar output is expected here!
            if any(fun_BC)
                try
                    bnd.val{ii}(fun_BC) = cellfun(@(f, coo) ...
                        {arrayfun(@(x, y) ...
                            f(x, y), coo(:, 1), coo(:, 2))}, ...
                        {bnd.val{ii}{fun_BC}}, ...
                        {bnd.bndDOF.bnd_DOF_coo{fun_BC}});
                catch ME
                   if strcmp(ME.identifier, 'MATLAB:arrayfun:NotAScalarOutput')
                       msg = ['Check bnd.val set up within the DRIVE '...
                           '- make shure that a function handle '...
                           'for the SCALAR FUNCTION is provided.'];
                       causeException = MException(...
                           'MATLAB:treatBC:wrongFunctionHandle', msg);
                       ME = addCause(ME, causeException);
                   end
                   rethrow(ME);
                end
            end

            % Replace basic info with bndDOF related values.
            if any(basic_BC)
                bnd.val{ii}(basic_BC) = cellfun(@(x, y) ...
                    {x + zeros(size(y))}, ...
                        {bnd.val{ii}{basic_BC}}, ...
                        {bnd.bndDOF.bnd_DOF{basic_BC}});
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

            case {'neumann', 'dtn'}
                % As Neumann (Robin) values needs to be evaluated at
                % quadrature nodes leave everything untouched.
                % See Fe.treatBC.m for its handling.
                if ~isfield(bnd, 'quad_ord')
                    warning(['field "quad_ord" - Quadrature order ', ...
                        'for bnd integral evaluation may be required. ', ...
                        'Set order to 1.']);
                    bnd.quad_ord = 1;
                end

            otherwise
                error('Bnd type: "%s" not supported yet.', bnd.type{ii});
        end
    end

    % Add parameter info (required in case of inhomogeneous Neumann BC).
    if any(strcmp(bnd.type, 'neumann'))
       bnd.param = param;
    end

    if verbosity
       fprintf('done.\n');
    end
end
