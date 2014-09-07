classdef Create_Functional_Map_Test < fcalign.test.Fc_Align_Test
% tests different usages of creating a functional map    

    methods (Test)
        
        function test_default(testCase)
        % tests default usage of creating a functional map
        
            % define parameters
            n = 16;
            p = 29;
            k = 9;
            dim = 6;
            t = 2;

            % make a positive definite matrix
            rng(0);
            data = rand(n, p);
            matrix = data * data';

            % create the expected map
            [V, l, d] = fcalign.internal.laplacian_eigs(matrix, k);
            expected = fcalign.create_diffusion_map(V, l, d, 'dim', dim, 't', t);
            expected.Ms_eig_vectors = V;
            expected.Ms_eig_values = l;
            
            % create the actual map
            actual = fcalign.create_functional_map(matrix, ...
                'number_eigs', k, 'map_dimension', dim, ...
                'diffusion_time', t);
            
            % verify that actual functional map is equal to expected one
            testCase.verify_equal_functional_map(actual, expected);
        end
        
    end
    
end
