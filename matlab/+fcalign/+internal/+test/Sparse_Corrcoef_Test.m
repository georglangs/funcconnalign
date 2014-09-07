classdef Sparse_Corrcoef_Test < fcalign.test.Fc_Align_Test
% test usage of creating a sparse correlation matrix    

    methods (Test)
        
        function test_default(obj)

            % define parameters
            n = 16;
            p = 29;
            threshold = 0.2;

            % reseed rng
            rng(0);

            % define some small random data
            data = rand(n, p);

            % create the expected sparse correlation matrix
            expected = corrcoef(data');
            expected = expected .* (expected > threshold);
            expected = sparse(expected);

            % create correlation matrix using method
            actual = fcalign.internal.sparse_corrcoef(data', threshold);

            % verify that actual is equal to expected
            obj.verifyEqual(actual, expected, 'RelTol', 0.01);
        end
        
    end
    
end
