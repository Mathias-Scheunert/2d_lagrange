classdef testRaviartThomas < matlab.unittest.TestCase
    % Tests Raviart-Thomas element implementation with MATLAB unittest.
    %
    % To run tests from project folder use:
    %   runtests('Test.testRaviartThomas')

    properties (TestParameter)
               
        factory = {[]};
    end

    methods (Test)
        
        function Basis(~)
            % Tests on basis (function) definitions.
            % TODO: implement. 
        end
        
        function Assembling(~)
            % Tests on assembling of mass and divergence matrix.
            % TODO: implement. 
        end
        
        function Mapping(~)
            % Tests on DOF mapping.
            % TODO: implement. 
        end
        
    end
end