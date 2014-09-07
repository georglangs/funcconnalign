classdef Create_Diffusion_Map_Test < fcalign.test.Fc_Align_Test
% tests different usages of creating a diffusion map    

    methods (Test)
        
        function test_default(testCase)
        % tests default usage of creating a functional map
        
            % define parameters
            n = 16;
            p = 29;
            k = 7;
            dim = 5;
            t = 1;
            
            % make a positive definite matrix
            rng(0);
            data = rand(n, p);
            matrix = data * data';
            
            % get the laplacian eigs
            [V, l, d] = fcalign.internal.laplacian_eigs(matrix, k);
            
            % create the expected sparse correlation matrix
            norm_factor = sqrt(sum(d));
            inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
            sqrt_D = spdiags(d.^0.5, 0, n, n);
            expected.Phi = sqrt_D*V(:, 1:dim) / norm_factor;
            expected.Psi = norm_factor * inv_sqrt_D*V(:, 1:dim);

            if expected.Psi(1, 1) < 0
                expected.Phi = expected.Phi * -1;
                expected.Psi = expected.Psi * -1;
            end
            
            Lt = spdiags(l(1:dim).^t, 0, dim, dim);
            expected.Gamma = expected.Psi*Lt;

            expected.diffusion_time = t;
            expected.map_dimension = dim;
            expected.delta = (l(dim) / l(2))^t;
            expected.determined_param = 'delta';
            
            % create correlation matrix using method
            actual = fcalign.create_diffusion_map(V, l, d, 'dim', dim, 't', t);

            % verify that expected correlation
            testCase.verify_equal_diffusion_map(actual, expected);
        end
        
    end
    
end
