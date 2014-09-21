classdef Align_Points_Test_Case < fcalign.test.Test_Case
% tests different usages of aligning multiple functional maps

    methods (Test)
        
        function test_default(this)
        % tests default usage of aligning multiple functional maps
        
            % define parameters
            n = 16;
            s = 32;
            
            % create a functional map
            maps(1) = this.create_functional_map('n', n, 's', s);
            
            % define a couple of transforms
            transform1 = fcalign.create_transform('seed', 1);
            transform2 = fcalign.create_transform('seed', 2);
            
            % apply these transforms to the map
            map2 = fcalign.transform_functional_map();
            map3 = fcalign.transform_functional_map();
            
            % align the maps
            alignedMaps
            
            % create the actual map
            actual = fcalign.create_functional_map(matrix, ...
                'number_eigs', k, 'map_dimension', dim, ...
                'diffusion_time', t);
            
            % verify that actual functional map is equal to expected one
            this.verify_equal_functional_map(actual, expected);
        end
        
    end
    
end
