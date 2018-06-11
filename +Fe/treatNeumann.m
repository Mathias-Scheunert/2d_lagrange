function [] = treatNewmann(mesh)
    % 
    %
    %
    % SYNTAX
    %  
    %
    % INPUT PARAMETER
    %   
    %
    % OUTPUT PARAMETER
    %   
    
    %% Check input.
    
    assert(isstruct(mesh) && all(isfield(mesh, {'cell2vtx', 'edge2vtx', 'vert2cord'})), ...
        'mesh - appended struct, including edge and mapping information, expected.');
    
    error('Not supported yet.');
    % TODO: implement.
end