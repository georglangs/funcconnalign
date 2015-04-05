function map = create_functional_map(W, varargin)
    % Creates a functional map from the given correlation matrix.
    %
    % Arguments
    % ---------
    % W : n x n numeric
    %   Correlation matrix.
    %
    % Keyword Arguments
    % -----------------
    % alpha : 0 <= float
    %   Value defining member of alpha family (default = 0).
    % n_eigs : 1 <= int <= n, t
    %   Number of laplacian eigenvalues to compute (default = n/2).
    % n_dims : 1 <= int <= n, t
    %   Number of dimension in functional map (default = 20)
    % center : bool
    %   Center laplacian matrix before taking eigendecomposition (default = false)
    % diff_time : 0 < numeric
    %   Diffusion time of functional map (default = 1)

    parser = inputParser();
    parser.addRequired('W', @(x) validateattributes(x, {'numeric'}, {'square'}));
    parser.addParamValue('alpha', 0, @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));
    parser.addParamValue('n_eigs', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    parser.addParamValue('n_dims', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    parser.addParamValue('center', false, @(x) validateattributes(x, {'logical'}));
    parser.addParamValue('diff_time', 1, @(x) validateattributes(x, {'numeric'}, {'positive'}));
    parser.addParamValue('verbose', false, @(x) validateattributes(x, {'logical'}));
    parser.parse(W, varargin{:});
    args = parser.Results;
    
    % process inputs further
    n = size(W, 1);
    if args.n_eigs == 0
        args.n_eigs = ceil(n * 0.5);
    end
    if args.n_dims == 0
        args.n_dims = ceil(args.n_eigs * 0.5);
    end
    
    % do some extra checks
    if args.n_eigs > n
        throw(fcalign.Exception('InvalidValue', [...
            'Number of eigenvalues (%u) is greater than size of ', ...
            'correlation matrix (%u)'], args.n_eigs, n));
    end

    if args.n_dims > args.n_eigs
        throw(fcalign.Exception('InvalidValue', [...
            'Number of map dimensions (%u) is greater than number of ', ...
            'eigenvalues (%u)'], args.n_dims, args.n_eigs));
    end

    % normalize correlation matrix by power of density according to family parameter
    if args.alpha > 0
        if args.verbose
            fprintf('Computing member of alpha family ...\n');
            tic_id = tic;
        end

        d = full(sum(W, 2));
        inv_D_alpha = spdiags(d.^-args.alpha, 0, n, n);
        W = inv_D_alpha*W*inv_D_alpha;

        if args.verbose
            fprintf('... done in %g s\n', toc(tic_id));
        end
    end

    % get eigenvectors of normalized laplacian
    if args.verbose
        fprintf('Taking eigendecomposition of normalized Laplacian ...\n');
        tic_id = tic;
    end

    [V, l, d] = fcalign.laplacian_eigs(W, args.n_eigs, 'center', args.center);

    if args.verbose
        fprintf('... done in %g s\n', toc(tic_id));
    end

    % create diffusion map
    map = fcalign.create_diffusion_map(V, l, d, 'dim', args.n_dims, 't', args.diff_time);
    map.Ms_eig_vectors = V;
    map.Ms_eig_values = l;
end
