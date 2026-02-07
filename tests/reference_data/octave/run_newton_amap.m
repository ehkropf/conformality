function [S, c, rho, a, history] = run_newton_amap(centers, radii, rotations, ...
                                                    c_ig, rho_ig, ...
                                                    N, newton_tol, max_newton_iter, ...
                                                    cgm_tol, max_cgm_iter, verbose)
% RUN_NEWTON_AMAP Run the full Newton iteration for the annulus case.
%   Replicates the newton_iteration loop from +bounded_map/bmap.m with amap methods.
%
%   [S, c, rho, a, history] = run_newton_amap(centers, radii, rotations,
%                                              c_ig, rho_ig,
%                                              N, newton_tol, max_newton_iter,
%                                              cgm_tol, max_cgm_iter, verbose)
%
%   Inputs:
%     centers, radii, rotations - target boundary parameters
%     c_ig, rho_ig              - initial guesses for canonical domain
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

    if nargin < 11, verbose = 1; end
    if nargin < 10, max_cgm_iter = 20; end
    if nargin < 9, cgm_tol = 1e-15; end
    if nargin < 8, max_newton_iter = 20; end
    if nargin < 7, newton_tol = 1e-14; end

    m = length(centers);

    % Initialize S (annulus sets c(1)=0)
    S = zeros(N, m);
    for j = 1:m
        S(:, j) = (0:N-1)' * 2*pi / N;
    end

    c = c_ig(:);
    c(1) = 0;  % annulus convention: first inner center is 0
    rho = rho_ig(:);

    history = struct('iteration', {}, 'normU', {}, 'D', {}, 'g', {}, ...
                     'U', {}, 'abs_eta', {}, 'S', {}, 'c', {}, 'rho', {});

    for iter = 1:max_newton_iter
        % Form system
        [D, g, abs_eta] = form_system_amap(S, c, rho, centers, radii, rotations, N);

        % Solve system
        U = solve_system_amap(D, g, N, cgm_tol, max_cgm_iter, verbose);

        % Newton update
        [S, c, rho, U_scaled] = newton_update_amap(S, c, rho, U, abs_eta, N);

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
            fprintf('amap iteration %d with update norm of %6e\n', iter, normU);
        end

        if normU <= newton_tol
            break;
        end
    end

    % Compute Fourier coefficients
    a = compute_fourier_coefficients(S, centers, radii, rotations, N);
end
