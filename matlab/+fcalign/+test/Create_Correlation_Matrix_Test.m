classdef Create_Correlation_Matrix_Test < fcalign.test.Fc_Align_Test
% test usage of creating a sparse correlation matrix    

    methods (Test)
        
        function test_default(testCase)
        % tests default usage of creating a correlation matrix
        
            % define parameters
            n = 16;
            t = 29;
            zero_inds = [3, 8];
            threshold = 0.2;

            % reseed rng
            rng(0);

            % define some small random fmri data
            data = rand(n, t);

            % make some of the rows be zeros
            data(zero_inds, :) = 0;

            % create the expected sparse correlation matrix
            expected.non_zeros = setdiff(1:n, [3,8])';
            non_zero_data = data(expected.non_zeros, :);
            expected.W = corrcoef(non_zero_data');
            expected.W = expected.W .* (expected.W > threshold);
            expected.W = sparse(expected.W);
            
            % create correlation matrix using method
            actual = fcalign.create_correlation_matrix(data, threshold);

            % verify that expected correlation
            testCase.verifyEqual(actual.W, expected.W, 'RelTol', 0.01);
            testCase.verifyEqual(actual.non_zeros, expected.non_zeros);
        end
        
    end
    
end
