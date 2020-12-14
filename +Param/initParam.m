function param = initParam(mesh, info)
    % Creates a parameter vector living on provided mesh.
    %
    % SYNTAX
    %   param = initParam(mesh, info)
    %
    % INPUT PARAMETER
    %   mesh ... Struct, containing the mesh information.
    %            For a detailed description of the content of the mesh
    %            struct please read header of Mesh.initMesh.
    %   info ... Struct, containing initial parameter name (Cell) and its
    %            conductivity value (Vector).
    %
    % OUTPUT PARAMETER
    %   param ... Vector [m x 1] of constant cell parameter values.

    %% Check input.

    assert(isstruct(mesh) && all(isfield(mesh, ...
            {'parameter_domain', 'parameter_domain_name'})), ...
        'mesh - Struct, including parameter domain mapping information, expected.');
    assert(isstruct(info) && all(isfield(info, ...
            {'val', 'name'})), ...
        'param - Struct, parameter domain names and corresponding values, expected.');

    %% Check consistency.

    assert(length(info.name) == length(info.val), ...
        'Mismatch between parameter domain names and value vector.');
    assert(length(mesh.parameter_domain_name) == length(info.name), ...
        ['Mismatch between number of parameter domains in mesh and ', ...
        'info initialized in DRIVE_. ', ...
        'Make sure, that all domains have physical names (Gmsh) and ', ...
        'that all domains are considered ind DRIVE_.']);
    if ~all(ismember(mesh.parameter_domain_name, info.name))
    names_err = {mesh.parameter_domain_name{:}; info.name{:}};
    text_err = [sprintf('mesh \t DRIVE_\n'), ...
        sprintf(repmat('%s \t %s\n', 1, length(names_err) / 2), ...
        names_err{:})];
    error(sprintf(['Unable to match physical domain names in mesh with ', ...
        'names provided in DRIVE_:\n', ...
        text_err]));
    end

    %% Map parameter domains to parameter vector.

    % Get assingement.
    [~, name_at_input] = ismember(mesh.parameter_domain_name, info.name);

    % Map parameter values to domain parts.
    domain2idx = mesh.parameter_domain == unique(mesh.parameter_domain).';
    param = domain2idx * info.val(name_at_input).';

end
