function [D, g, abs_eta] = form_system_ellipse(S, c, rho, norm_cond, ...
                                                bdry_centers, bdry_semi_maj, ...
                                                bdry_semi_min, bdry_rotation, N)
% FORM_SYSTEM_ELLIPSE Form the linear system for the non-annulus case with ellipse boundaries.
%   Like form_system_namap but uses ellipse_xi_eta instead of circular_xi_eta.
%
%   [D, g, abs_eta] = form_system_ellipse(S, c, rho, norm_cond,
%                                          bdry_centers, bdry_semi_maj,
%                                          bdry_semi_min, bdry_rotation, N)
%
%   Inputs:
%     S              - boundary parameters, N x m matrix
%     c              - circle centers in canonical domain, (m-1) x 1
%     rho            - circle radii in canonical domain, (m-1) x 1
%     norm_cond      - normalization conditions [alpha, z0, f0], 3 x 1
%     bdry_centers   - target boundary centers, m x 1 (complex)
%     bdry_semi_maj  - target boundary semi-major axes, m x 1
%     bdry_semi_min  - target boundary semi-minor axes, m x 1
%     bdry_rotation  - target boundary rotation angles, m x 1
%     N              - number of boundary points per boundary
%
%   Outputs:
%     D       - system matrix
%     g       - right-hand side vector
%     abs_eta - absolute value of tangent vectors, N x m

    m = size(S, 2);
    M = N/2;

    q = exp(-1i*(0:N-1)'*2*pi/N);
    tl = 2*pi * ones(1, m);  % bellipse.tot_len = 2*pi

    % Compute xi, eta for all boundaries using ellipse parameterization
    xi = zeros(N, m);
    eta = zeros(N, m);
    for j = 1:m
        [xi(:,j), eta(:,j)] = ellipse_xi_eta(bdry_centers(j), ...
                                               bdry_semi_maj(j), ...
                                               bdry_semi_min(j), ...
                                               bdry_rotation(j), S(:,j));
    end

    abs_eta = abs(eta);
    eta = eta ./ abs_eta;

    % Size of U vector
    U_row_sz = m*N + 3*(m-1);
    D = zeros(m*M + 2, U_row_sz);

    % p1 row vector
    p1 = [1 zeros(1, N-1)];
    if norm_cond(2) ~= 0
        p1(1, 2:M) = norm_cond(2).^(1:M-1);
    end

    % nu = 1 (outer boundary)
    Pnu = [make_Pnu(1, m, c, rho, N); p1];
    D(1:end-1, 1:N) = Pnu * fft(diag(eta(:,1)));
    g = -Pnu * fft(xi(:,1));

    % pnu row vector (reused for each inner boundary)
    pnu = zeros(1, N);

    for nu = 2:m
        pnu(M+1:N) = (rho(nu-1) / (norm_cond(2) - c(nu-1))).^(M:-1:1);
        Pnu = [make_Pnu(nu, m, c, rho, N); pnu];
        D(1:end-1, (nu-1)*N+1:nu*N) = Pnu * fft(diag(eta(:,nu)));

        % W columns
        Sdiff = diff([S(:,nu); S(1,nu) + tl(nu)]) * (N/(2*pi));
        zeta = 1i * (abs_eta(:,nu) .* eta(:,nu) .* Sdiff) / rho(nu-1);
        D(1:end-1, m*N+nu-1) = Pnu * fft(zeta);
        D(1:end-1, m*N+m+2*(nu-2)) = Pnu * fft(q .* zeta);
        D(1:end-1, m*N+m+2*(nu-2)+1) = 1i * D(1:end-1, m*N+m+2*(nu-2));

        g = g - Pnu * fft(xi(:,nu));
    end

    % Add normalization conditions
    D(end, :) = [1 zeros(1, m*N + 3*(m-1) - 1)];
    g(end) = g(end) + N * norm_cond(3);
    g = [g; 0];
end
