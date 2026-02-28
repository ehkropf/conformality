% THESIS3_ITERATION_DEBUG Generate detailed per-iteration reference data
% for the thesis example 3 (identity map, m=4, N=256) to diagnose
% Newton iteration divergence in the C++ implementation.
%
% Generates JSON files capturing full state at each stage across two
% Newton iterations:
%   - Initial state (S, c, rho)
%   - After formSystem (D, g, abs_eta) for iterations 1 and 2
%   - After solveSystem (U) for iterations 1 and 2
%   - After newtonUpdate (S, c, rho, normU) for iterations 1 and 2
%
% Run: octave --no-gui tests/reference_data/octave/thesis3_iteration_debug.m

output_dir = 'tests/reference_data/data';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

fprintf('=== Thesis3 iteration debug: m=4, N=256 ===\n');

N = 256;
centers = [0; -0.5; 0.25+0.43i; 0.25-0.43i];
radii = [1; 0.25; 0.25; 0.25];
rotations = [0; 0; 0; 0];
c_ig = [-0.4; 0.35+0.43i; 0.35-0.43i];
rho_ig = [0.25; 0.25; 0.25];
norm_cond = [1; 0; 0];
m = 4;

cgm_tol = 1e-15;
max_cgm_iter = 100;
verbose = 1;

% Initialize S
S = zeros(N, m);
for j = 1:m
    S(:, j) = (0:N-1)' * 2*pi / N;
end
c = c_ig(:);
rho = rho_ig(:);

% --- Save initial state ---
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'initial');
data.metadata.centers_re = real(centers).';
data.metadata.centers_im = imag(centers).';
data.metadata.radii = radii.';
data.metadata.rotations = rotations.';
data.metadata.norm_cond = norm_cond.';
data.matrices = struct('S', S);
data.vectors = struct('c', c, 'rho', rho);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'thesis3_debug_initial.json'), data);
fprintf('  Wrote thesis3_debug_initial.json\n');

% ==================== Iteration 1 ====================
fprintf('\n--- Iteration 1 ---\n');

% Form system
[D1, g1, abs_eta1] = form_system_namap(S, c, rho, norm_cond, ...
                                        centers, radii, rotations, N);

% Save form_system results for iteration 1
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'form_system', 'iteration', 1);
data.matrices = struct('D', D1, 'abs_eta', abs_eta1);
data.vectors = struct('g', g1);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'thesis3_debug_form_system_iter1.json'), data);
fprintf('  Wrote thesis3_debug_form_system_iter1.json\n');

% Solve system
U1 = solve_system_namap(D1, g1, N, cgm_tol, max_cgm_iter, verbose);

% Save solve_system results for iteration 1
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'solve_system', 'iteration', 1);
data.vectors = struct('U', U1);
data.scalars = struct();
data.matrices = struct();
save_reference_json(fullfile(output_dir, 'thesis3_debug_solve_system_iter1.json'), data);
fprintf('  Wrote thesis3_debug_solve_system_iter1.json\n');

% Newton update
[S, c, rho, U1_scaled] = newton_update_namap(S, c, rho, U1, abs_eta1, N);
normU1 = norm(U1_scaled, inf);
fprintf('  Iteration 1 update norm: %e\n', normU1);

% Save newton_update results for iteration 1
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'newton_update', 'iteration', 1);
data.matrices = struct('S', S, 'abs_eta', abs_eta1);
data.vectors = struct('c', c, 'rho', rho, 'U_scaled', U1_scaled);
data.scalars = struct('normU', normU1);
save_reference_json(fullfile(output_dir, 'thesis3_debug_newton_update_iter1.json'), data);
fprintf('  Wrote thesis3_debug_newton_update_iter1.json\n');

% ==================== Iteration 2 ====================
fprintf('\n--- Iteration 2 ---\n');

% Form system with updated state
[D2, g2, abs_eta2] = form_system_namap(S, c, rho, norm_cond, ...
                                        centers, radii, rotations, N);

% Check for inf/nan in D2 and g2
fprintf('  D2 max abs: %e\n', max(abs(D2(:))));
fprintf('  g2 max abs: %e\n', max(abs(g2(:))));
fprintf('  D2 has inf: %d, has nan: %d\n', any(isinf(D2(:))), any(isnan(D2(:))));
fprintf('  g2 has inf: %d, has nan: %d\n', any(isinf(g2(:))), any(isnan(g2(:))));

% Save form_system results for iteration 2
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'form_system', 'iteration', 2);
data.matrices = struct('D', D2, 'abs_eta', abs_eta2);
data.vectors = struct('g', g2);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'thesis3_debug_form_system_iter2.json'), data);
fprintf('  Wrote thesis3_debug_form_system_iter2.json\n');

% Solve system (iteration 2)
U2 = solve_system_namap(D2, g2, N, cgm_tol, max_cgm_iter, verbose);

% Save solve_system results for iteration 2
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'solve_system', 'iteration', 2);
data.vectors = struct('U', U2);
data.scalars = struct();
data.matrices = struct();
save_reference_json(fullfile(output_dir, 'thesis3_debug_solve_system_iter2.json'), data);
fprintf('  Wrote thesis3_debug_solve_system_iter2.json\n');

% Newton update (iteration 2)
[S, c, rho, U2_scaled] = newton_update_namap(S, c, rho, U2, abs_eta2, N);
normU2 = norm(U2_scaled, inf);
fprintf('  Iteration 2 update norm: %e\n', normU2);

% Save newton_update results for iteration 2
data = struct();
data.metadata = struct('case', 'thesis3_debug', 'N', N, 'm', m, ...
                       'stage', 'newton_update', 'iteration', 2);
data.matrices = struct('S', S, 'abs_eta', abs_eta2);
data.vectors = struct('c', c, 'rho', rho, 'U_scaled', U2_scaled);
data.scalars = struct('normU', normU2);
save_reference_json(fullfile(output_dir, 'thesis3_debug_newton_update_iter2.json'), data);
fprintf('  Wrote thesis3_debug_newton_update_iter2.json\n');

fprintf('\n=== Thesis3 debug data generation complete ===\n');
