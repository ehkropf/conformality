function [S, c, rho, a, history] = run_newton_namap(centers, radii, rotations, ...
                                                     c_ig, rho_ig, norm_cond, ...
                                                     N, newton_tol, max_newton_iter, ...
                                                     cgm_tol, max_cgm_iter, verbose)
% RUN_NEWTON_NAMAP Run the full Newton iteration for the non-annulus case.
%   Replicates the newton_iteration loop from +bounded_map/bmap.m with namap methods.
%
%   [S, c, rho, a, history] = run_newton_namap(centers, radii, rotations,
%                                               c_ig, rho_ig, norm_cond,
%                                               N, newton_tol, max_newton_iter,
%                                               cgm_tol, max_cgm_iter, verbose)
%
%   Inputs:
%     centers, radii, rotations - target boundary parameters
%     c_ig, rho_ig              - initial guesses for canonical domain
%     norm_cond                 - normalization conditions [alpha, z0, f0]
%     N                         - number of boundary points per boundary
%     newton_tol                - Newton convergence tolerance
%     max_newton_iter           - maximum Newton iterations
%     cgm_tol                   - CG solver tolerance
%     max_cgm_iter              - maximum CG iterations
%     verbose                   - verbosity flag
%
%   Outputs:
%     S       - final boundary parameters
%     c       - final canonical domain centers
%     rho     - final canonical domain radii
%     a       - Fourier coefficients
%     history - struct array with per-iteration data

    if nargin < 12, verbose = 1; end
    if nargin < 11, max_cgm_iter = 20; end
    if nargin < 10, cgm_tol = 1e-15; end
    if nargin < 9, max_newton_iter = 20; end
    if nargin < 8, newton_tol = 1e-14; end

    m = length(centers);

    % Initialize S
    S = zeros(N, m);
    for j = 1:m
        S(:, j) = (0:N-1)' * 2*pi / N;
    end

    c = c_ig(:);
    rho = rho_ig(:);

    history = struct('iteration', {}, 'normU', {}, 'D', {}, 'g', {}, ...
                     'U', {}, 'abs_eta', {}, 'S', {}, 'c', {}, 'rho', {});

    for iter = 1:max_newton_iter
        % Form system
        [D, g, abs_eta] = form_system_namap(S, c, rho, norm_cond, ...
                                             centers, radii, rotations, N);

        % Solve system
        U = solve_system_namap(D, g, N, cgm_tol, max_cgm_iter, verbose);

        % Newton update
        [S, c, rho, U_scaled] = newton_update_namap(S, c, rho, U, abs_eta, N);

        normU = norm(U_scaled, inf);

        % Save history
        h = struct();
        h.iteration = iter;
        h.normU = normU;
        h.D = D;
        h.g = g;
        h.U = U;
        h.abs_eta = abs_eta;
        h.S = S;
        h.c = c;
        h.rho = rho;
        history(iter) = h;

        if verbose
            fprintf('namap iteration %d with update norm of %6e\n', iter, normU);
        end

        if normU <= newton_tol
            break;
        end
    end

    % Compute Fourier coefficients
    a = compute_fourier_coefficients(S, centers, radii, rotations, N);
end
