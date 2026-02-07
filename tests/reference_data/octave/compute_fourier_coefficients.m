function a = compute_fourier_coefficients(S, centers, radii, rotations, N)
% COMPUTE_FOURIER_COEFFICIENTS Compute Fourier coefficients of the map.
%   Replicates set_fcoef from +bounded_map/namap.m and amap.m
%
%   a = compute_fourier_coefficients(S, centers, radii, rotations, N)
%
%   Inputs:
%     S         - boundary parameters, N x m matrix
%     centers   - target boundary centers, m x 1
%     radii     - target boundary radii, m x 1
%     rotations - target boundary rotations, m x 1
%     N         - number of boundary points per boundary
%
%   Outputs:
%     a - Fourier coefficients, N x m matrix

    m = size(S, 2);

    xi = zeros(N, m);
    for nu = 1:m
        xi(:, nu) = circular_xi_eta(centers(nu), radii(nu), rotations(nu), ...
                                     S(:, nu));
    end

    a = zeros(N, m);
    for nu = 1:m
        a(:, nu) = fft(xi(:, nu)) / N;
    end
end
