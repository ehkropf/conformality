function [D, g, abs_eta] = form_system_amap(S, c, rho, ...
                                            centers, radii, rotations, N)
% FORM_SYSTEM_AMAP Form the linear system for the annulus case.
%   Replicates amap.form_system from +bounded_map/amap.m
%
%   [D, g, abs_eta] = form_system_amap(S, c, rho,
%                                       centers, radii, rotations, N)
%
%   Inputs:
%     S          - boundary parameters, N x m matrix
%     c          - circle centers in canonical domain, (m-1) x 1
%     rho        - circle radii in canonical domain, (m-1) x 1
%     centers    - target boundary centers, m x 1
%     radii      - target boundary radii, m x 1
%     rotations  - target boundary rotations, m x 1
%     N          - number of boundary points per boundary
%
%   Outputs:
%     D       - system matrix
%     g       - right-hand side vector
%     abs_eta - absolute value of tangent vectors, N x m

    m = size(S, 2);
    M = N/2;

    q = exp(-1i*(0:N-1)'*2*pi/N);
    tl = 2*pi * ones(1, m);  % total length for circles is always 2*pi

    % Compute xi, eta for all boundaries
    xi = zeros(N, m);
    eta = zeros(N, m);
    for j = 1:m
        [xi(:,j), eta(:,j)] = circular_xi_eta(centers(j), radii(j), ...
                                                rotations(j), S(:,j));
    end

    abs_eta = abs(eta);
    eta = eta ./ abs_eta;

    % Size of U vector
    U_row_sz = m*N + 1 + (m >= 3)*2 + (m >= 4)*3*(m-3);
    D = zeros(m*M, U_row_sz);

    % nu = 1 (outer boundary)
    Pnu = make_Pnu_ann(1, m, c, rho, N);
    D(:, 1:N) = Pnu * fft(diag(eta(:,1)));
    g = -Pnu * fft(xi(:,1));

    for nu = 2:m
        Pnu = make_Pnu_ann(nu, m, c, rho, N);
        D(:, (nu-1)*N+1:nu*N) = Pnu * fft(diag(eta(:,nu)));

        Sdiff = diff([S(:,nu); S(1,nu) + tl(nu)]);
        zeta = 1i * (abs_eta(:,nu) .* eta(:,nu) .* (Sdiff * (N/(2*pi)))) ...
               / rho(nu-1);
        D(:, m*N+nu-1) = Pnu * fft(zeta);
        if nu == 3
            D(:, m*N+m) = Pnu * fft(q .* zeta);
        elseif nu > 3
            D(:, m*N+m+2*(nu-3)-1) = Pnu * fft(q .* zeta);
            D(:, m*N+m+2*(nu-3)) = 1i * D(:, m*N+m+2*(nu-3)-1);
        end

        g = g - Pnu * fft(xi(:,nu));
    end
end
