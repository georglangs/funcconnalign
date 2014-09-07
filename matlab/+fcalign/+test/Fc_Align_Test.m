classdef (Abstract) Fc_Align_Test < matlab.unittest.TestCase
% test usage of a component of the fc align package

    methods (Test, Abstract)
       
        test_default(obj);
        % tests default usage
       
    end
    
    methods (TestMethodSetup)
        
        function test_method_setup(obj)
        % setup for each test method
            rng(0);
        end
        
    end
            
    
    methods (Access = protected)
        
        function verify_equal_diffusion_map(obj, actual, expected)
        % verifies that the actual diffusion map is equal to the expected one
            obj.assertClass(actual, 'struct');
            obj.assertClass(expected, 'struct');
        
            obj.verifyEqual(actual.Phi, expected.Phi, 'RelTol', 0.01);
            obj.verifyEqual(actual.Psi, expected.Psi, 'RelTol', 0.01);
            obj.verifyEqual(actual.Gamma, expected.Gamma, 'RelTol', 0.01);
            obj.verifyEqual(actual.diffusion_time, expected.diffusion_time, 'RelTol', 0.01);
            obj.verifyEqual(actual.map_dimension, expected.map_dimension, 'RelTol', 0.01);
            obj.verifyEqual(actual.delta, expected.delta, 'RelTol', 0.01);
            obj.verifyEqual(actual.determined_param, expected.determined_param, 'RelTol', 0.01);
        end
        
        function verify_equal_functional_map(obj, actual, expected)
        % verifies that the actual functional map is equal to the expected one
            obj.assertClass(actual, 'struct');
            obj.assertClass(expected, 'struct');
        
            obj.verifyEqual(actual.Ms_eig_vectors, expected.Ms_eig_vectors, 'RelTol', 0.01);
            obj.verifyEqual(actual.Ms_eig_values, expected.Ms_eig_values, 'RelTol', 0.01);
            obj.verify_equal_diffusion_map(actual, expected);
        end
        
    end
    
end
