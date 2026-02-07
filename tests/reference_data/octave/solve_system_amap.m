function U = solve_system_amap(D, g, N, cgm_tol, max_cgm_iter, verbose)
% SOLVE_SYSTEM_AMAP Solve the normal equations via CG for annulus case.
%   Replicates amap.solve_system from +bounded_map/amap.m
%
%   U = solve_system_amap(D, g, N, cgm_tol, max_cgm_iter, verbose)
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

    % Annulus case forms explicit matrix (not function handle)
    A = 2*real(D'*D)/N;
    b = 2*real(D'*g)/N;

    U = cgm(A, b, cgm_tol, max_cgm_iter, verbose);
end
