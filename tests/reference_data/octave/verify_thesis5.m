% VERIFY_THESIS5 Run Fornberg MC on thesis example 5 and output reference moduli.
%
% Thesis example 5: m=3, ellipse outer boundary, two inner ellipses.
% Domain setup from design/fornberg/fornmc/th_gen_ex5.m
%
% Run from the project root:
%   octave --no-gui --path design/fornberg/fornmc \
%     tests/reference_data/octave/verify_thesis5.m
%
% Or with explicit path:
%   octave --no-gui --path /path/to/design/fornberg/fornmc \
%     tests/reference_data/octave/verify_thesis5.m
%
% Requires make_Pnu.m from design/fornberg/fornmc on the path.
% Outputs JSON to tests/reference_data/data/thesis5_moduli.json

% Add the octave scripts directory to the path
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

% Try to find design/fornberg/fornmc for make_Pnu
design_dir = fullfile(script_dir, '..', '..', '..', 'design', 'fornberg', 'fornmc');
if exist(design_dir, 'dir')
    addpath(design_dir);
end

% Verify make_Pnu is available
if ~exist('make_Pnu', 'file')
    error(['make_Pnu.m not found. Add design/fornberg/fornmc to path:\n' ...
           '  octave --no-gui --path design/fornberg/fornmc %s'], ...
          mfilename('fullpath'));
end

output_dir = fullfile(script_dir, '..', 'data');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% ==================== Thesis Example 5 Setup ====================
% From th_gen_ex5.m:
%   C1 = bellipse([0    2     3/2],       N)
%   C2 = bellipse([-.8  3/16  3/8   0],   N, [-.47  .15])
%   C3 = bellipse([.7-.1i  3/8  3/16  pi/4], N, [.4-.1i  .15])

N = 2^8;  % 256, same as thesis
m = 3;

% Boundary parameters: each row is [center, semi_maj, semi_min, rotation]
bdry_centers  = [0; -0.8; 0.7-0.1i];
bdry_semi_maj = [2; 3/16; 3/8];
bdry_semi_min = [3/2; 3/8; 3/16];
bdry_rotation = [0; 0; pi/4];

% Initial guesses for canonical domain (inner boundaries only)
c_ig  = [-0.47; 0.4-0.1i];
rho_ig = [0.15; 0.15];

% Normalization conditions [alpha, z0, f0]
norm_cond = [1; 0; 0];

% Algorithm parameters
newton_tol = 1e-14;
max_newton_iter = 50;
cgm_tol = 1e-15;
max_cgm_iter = 20;

fprintf('\n=== Thesis Example 5: m=%d, N=%d ===\n', m, N);
fprintf('Boundaries:\n');
for j = 1:m
    fprintf('  C%d: center=%g%+gi, a=%g, b=%g, rot=%g\n', ...
            j, real(bdry_centers(j)), imag(bdry_centers(j)), ...
            bdry_semi_maj(j), bdry_semi_min(j), bdry_rotation(j));
end
fprintf('Initial guesses:\n');
for j = 1:m-1
    fprintf('  c_%d = %g%+gi, rho_%d = %g\n', ...
            j+1, real(c_ig(j)), imag(c_ig(j)), j+1, rho_ig(j));
end

%% ==================== Newton Iteration ====================
% Initialize S (uniform parameter distribution)
S = zeros(N, m);
for j = 1:m
    S(:, j) = (0:N-1)' * 2*pi / N;
end

c = c_ig(:);
rho = rho_ig(:);

for iter = 1:max_newton_iter
    % Form system
    [D, g, abs_eta] = form_system_ellipse(S, c, rho, norm_cond, ...
                                           bdry_centers, bdry_semi_maj, ...
                                           bdry_semi_min, bdry_rotation, N);

    % Solve system (CG on normal equations)
    U = solve_system_namap(D, g, N, cgm_tol, max_cgm_iter, 0);

    % Newton update
    [S, c, rho, U_scaled] = newton_update_namap(S, c, rho, U, abs_eta, N);

    normU = norm(U_scaled, inf);
    fprintf('  iteration %d: ||U|| = %e\n', iter, normU);

    if normU <= newton_tol
        fprintf('Converged at iteration %d.\n', iter);
        break;
    end
end

num_iterations = iter;

%% ==================== Display Results ====================
fprintf('\n=== Conformal Moduli ===\n');
for j = 1:m-1
    fprintf('  c_%d = (%.15g, %.15g)  rho_%d = %.15g\n', ...
            j+1, real(c(j)), imag(c(j)), j+1, rho(j));
end

%% ==================== C++ Comparison ====================
fprintf('\n=== C++ Results (from GH-106) ===\n');
fprintf('  c_2 = (-0.470892, -0.00839284)  rho_2 = 0.150536\n');
fprintf('  c_3 = (0.399392, -0.0862608)    rho_3 = 0.155272\n');

fprintf('\n=== Differences ===\n');
cpp_c = [-0.470892 - 0.00839284i; 0.399392 - 0.0862608i];
cpp_rho = [0.150536; 0.155272];
for j = 1:m-1
    dc = abs(c(j) - cpp_c(j));
    dr = abs(rho(j) - cpp_rho(j));
    fprintf('  |delta c_%d| = %e  |delta rho_%d| = %e\n', j+1, dc, j+1, dr);
end

%% ==================== Save JSON ====================
% Compute Fourier coefficients
a = zeros(N, m);
for nu = 1:m
    xi_nu = ellipse_xi_eta(bdry_centers(nu), bdry_semi_maj(nu), ...
                            bdry_semi_min(nu), bdry_rotation(nu), S(:,nu));
    a(:, nu) = fft(xi_nu) / N;
end

data = struct();
data.metadata = struct('case', 'thesis5_ellipse_m3', ...
                       'N', N, ...
                       'm', m, ...
                       'stage', 'converged', ...
                       'num_iterations', num_iterations, ...
                       'newton_tol', newton_tol);
data.metadata.bdry_centers_re = real(bdry_centers).';
data.metadata.bdry_centers_im = imag(bdry_centers).';
data.metadata.bdry_semi_maj = bdry_semi_maj.';
data.metadata.bdry_semi_min = bdry_semi_min.';
data.metadata.bdry_rotation = bdry_rotation.';
data.matrices = struct('S', S, 'a', a);
data.vectors = struct('c', c, 'rho', rho);
data.scalars = struct('final_normU', normU);

output_file = fullfile(output_dir, 'thesis5_moduli.json');
save_reference_json(output_file, data);
fprintf('\nWrote %s\n', output_file);

fprintf('\n=== Done ===\n');
