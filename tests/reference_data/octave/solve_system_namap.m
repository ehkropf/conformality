function U = solve_system_namap(D, g, N, cgm_tol, max_cgm_iter, verbose)
% SOLVE_SYSTEM_NAMAP Solve the normal equations via CG for non-annulus case.
%   Replicates namap.solve_system from +bounded_map/namap.m
%
%   U = solve_system_namap(D, g, N, cgm_tol, max_cgm_iter, verbose)
%
%   Inputs:
%     D            - system matrix
%     g            - right-hand side vector
%     N            - number of boundary points per boundary
%     cgm_tol      - CG solver tolerance
%     max_cgm_iter - maximum CG iterations
%     verbose      - verbosity flag
%
%   Outputs:
%     U - solution vector

    if nargin < 6, verbose = 1; end
    if nargin < 5, max_cgm_iter = 20; end
    if nargin < 4, cgm_tol = 1e-15; end

    b = 2*real(D'*g)/N;

    % Use function handle for matrix-vector product (matches MATLAB namap)
    Amx = @(x) 2*(real(D)'*(real(D)*x) + imag(D)'*(imag(D)*x))/N;

    U = cgm(Amx, b, cgm_tol, max_cgm_iter, verbose);
end
