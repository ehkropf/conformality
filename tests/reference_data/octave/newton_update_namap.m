function [S_new, c_new, rho_new, U_scaled] = newton_update_namap(S, c, rho, ...
                                                                  U, abs_eta, N)
% NEWTON_UPDATE_NAMAP Apply Newton update for non-annulus case.
%   Replicates namap.newton_update from +bounded_map/namap.m
%
%   [S_new, c_new, rho_new, U_scaled] = newton_update_namap(S, c, rho,
%                                                            U, abs_eta, N)
%
%   Inputs:
%     S        - boundary parameters, N x m matrix
%     c        - circle centers in canonical domain, (m-1) x 1
%     rho      - circle radii in canonical domain, (m-1) x 1
%     U        - solution vector from solve_system
%     abs_eta  - absolute value of tangent vectors, N x m
%     N        - number of boundary points per boundary
%
%   Outputs:
%     S_new    - updated boundary parameters
%     c_new    - updated circle centers
%     rho_new  - updated circle radii
%     U_scaled - U after scaling by abs_eta (for convergence check)

    m = size(S, 2);

    % Scale S-update by 1/abs_eta
    U_scaled = [U(1:m*N) ./ abs_eta(:); U(m*N+1:end)];

    S_new = S + reshape(U_scaled(1:m*N), N, m);
    rho_new = rho + U_scaled(m*N+1:m*N+m-1);
    c_new = c + U_scaled(m*N+m:2:m*N+m+2*m-4) ...
              + 1i*U_scaled(m*N+m+1:2:m*N+m+2*m-2);
end
