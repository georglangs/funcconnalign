classdef (Abstract) Test_Case < matlab.unittest.TestCase
% test usage of a component of the fc align package
    
    methods (Access = protected)
        
        function data = create_data(this, varargin)
        % creates a standard data matrix
        %
        % OPTIONAL INPUTS
        % n             number of nodes
        %               1 <= integer (default = 16)
        % s             number of samples
        %               1 <= integer (default = 32)
        % seed          seed for random number generator
        %               0 <= integer (default = 0)
        
            % parse inputs
            parser = inputParser();
            parser.addParamValue('n', 16, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('s', 32, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.parse(varargin{:});
            inputs = parser.Results;
        
            % generate data from random uniform
            rng(inputs.seed);
            data = rand(inputs.n, inputs.s);
        end
        
        function matrix = create_correlation_matrix(this, varargin)
        % creates a standard correlation matrix
        %
        % OPTIONAL INPUTS
        % n             number of nodes
        %               1 <= integer (default = 16)
        % s             number of samples
        %               1 <= integer (default = 32)
        % threshold     correlation threshold
        %               -1 <= double <= 1 (default = 0.1)
        % seed          seed for random number generator
        %               0 <= integer (default = 0)
        
            % parse inputs
            parser = inputParser();
            parser.addParamValue('n', 16, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('s', 32, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('threshold', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -1, '<=', 1}));
            parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.parse(varargin{:});
            inputs = parser.Results;
        
            % generate data then calculate correlation matrix
            data = this.create_data('n', inputs.n, 's', inputs.s, 'seed', inputs.seed);
            matrix = fcalign.create_correlation_matrix(data, inputs.threshold);
        end
        
        function map = create_functional_map(this, varargin)
        % creates a standard functional map
        %
        % OPTIONAL INPUTS
        % n             number of nodes
        %               1 <= integer (default = 16)
        % s             number of samples
        %               1 <= integer (default = 32)
        % threshold     correlation threshold
        %               -1 <= double <= 1 (default = 0.1)
        % seed          seed for random number generator
        %               0 <= integer (default = 0)
        
            % parse inputs
            parser = inputParser();
            parser.addParamValue('n', 16, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('s', 32, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('threshold', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -1, '<=', 1}));
            parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.parse(varargin{:});
            inputs = parser.Results;
                
            % create map from correlation matrix
            matrix = this.create_correlation_matrix('n', inputs.n, 's', inputs.s, 'threshold', inputs.threshold, 'seed', inputs.seed);
            map = fcalign.create_functional_map(matrix);
        end
        
        function transform = create_transform(this, varargin)
        % creates a p-dimension transform
        %
        % OPTIONAL INPUTS
        % p         dimension of transform
        %           1 <= integer (default = 20)
        % type      type of transformation
        %           string from (default = rigid)
        %               - 'rotation': determinant is 1
        %               - 'orthonormal': allow flips (determinant can be -1 or +1)
        %               - 'rigid': rotation and translation
        %               - 'scale_orthonormal': incorporates global scale factor
        % seed      seed for random number generator
        %           0 <= integer (default = 0)
        
            % parse inputs
            parser = inputParser();
            parser.addParamValue('p', 20, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.addParamValue('type', 'rigid', @(x) validatestring(x, {'rotation', 'orthonormal', 'rigid', 'scale'}));
            parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
            parser.parse(varargin{:});
            inputs = parser.Results;
                
            % create a skew symmetric matrix
            skewSymmetric = rand(p, p);
            for i = 1:p
                for j = (i+1):p
                    skewSymmetric(i, j) = -skewSymmetric(i, j);
                end
            end
            
            % exponentiate it to get a rotation
            rotationMatrix = expm(skewSymmetric);
            
            % create random 
            rng(inputs.seed);
            switch inputs.type
                case 'rotation'
                    transform = rotationMatrix;
                case 'orthonormal'
                    [V, L] = eig(rotationMatrix);
                    L(1,1) = -L(1, 1);
                    transform = V*L*V';
                case 'rigid'
                    transform = ;
                case 'scale'
                    transform = ;
                otherwise
                    fcalign.internal.throwException('', '');
            end
                
        end
        
        function verify_equal_diffusion_map(this, actual, expected)
        % verifies that the actual diffusion map is equal to the expected one
        %
        % REQUIRED INPUTS
        % actual        actual diffusion map
        %               struct
        % expected      expected diffusion map
        %               struct
            this.assertClass(actual, 'struct');
            this.assertClass(expected, 'struct');
        
            this.verifyEqual(actual.Phi, expected.Phi, 'RelTol', 0.01);
            this.verifyEqual(actual.Psi, expected.Psi, 'RelTol', 0.01);
            this.verifyEqual(actual.Gamma, expected.Gamma, 'RelTol', 0.01);
            this.verifyEqual(actual.diffusion_time, expected.diffusion_time, 'RelTol', 0.01);
            this.verifyEqual(actual.map_dimension, expected.map_dimension, 'RelTol', 0.01);
            this.verifyEqual(actual.delta, expected.delta, 'RelTol', 0.01);
            this.verifyEqual(actual.determined_param, expected.determined_param, 'RelTol', 0.01);
        end
        
        function verify_equal_functional_map(this, actual, expected)
        % verifies that the actual functional map is equal to the expected one
        %
        % REQUIRED INPUTS
        % actual        actual functional map
        %               struct
        % expected      expected functional map
        %               struct
            this.assertClass(actual, 'struct');
            this.assertClass(expected, 'struct');
        
            this.verifyEqual(actual.Ms_eig_vectors, expected.Ms_eig_vectors, 'RelTol', 0.01);
            this.verifyEqual(actual.Ms_eig_values, expected.Ms_eig_values, 'RelTol', 0.01);
            this.verify_equal_diffusion_map(actual, expected);
        end
        
    end
    
end
