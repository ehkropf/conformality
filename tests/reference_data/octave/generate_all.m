% GENERATE_ALL Generate all MATLAB reference data for C++ comparison tests.
%
% Run from the project root:
%   octave --no-gui --path tests/reference_data/octave --path design/fornberg/fornmc \
%     --eval "generate_all"
%
% Generates JSON files in tests/reference_data/data/

output_dir = 'tests/reference_data/data';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

newton_tol = 1e-14;
max_newton_iter = 20;
cgm_tol = 1e-15;
max_cgm_iter = 20;
verbose = 0;

%% ==================== Case 1: Annulus (m=2, N=64) ====================
fprintf('=== Generating annulus (m=2) reference data ===\n');

N = 64;
centers = [0; 0.3];
radii = [1.0; 0.15];
rotations = [0; 0];
c_ig = [0.3];
rho_ig = [0.15];
m = 2;

% --- P matrices ---
S0 = zeros(N, m);
for j = 1:m
    S0(:, j) = (0:N-1)' * 2*pi / N;
end

P_ann = cell(m, 1);
for nu = 1:m
    P_ann{nu} = make_Pnu_ann(nu, m, c_ig, rho_ig, N);
end

for nu = 1:m
    data = struct();
    data.metadata = struct('case', 'annulus', 'N', N, 'm', m, 'nu', nu, ...
                           'stage', 'P_matrix');
    data.matrices = struct();
    data.matrices.P = P_ann{nu};
    data.vectors = struct();
    data.scalars = struct();
    save_reference_json(fullfile(output_dir, ...
        sprintf('annulus_P_nu%d.json', nu)), data);
    fprintf('  Wrote annulus_P_nu%d.json\n', nu);
end

% --- Full Newton iteration with history ---
[S_final, c_final, rho_final, a_final, history] = ...
    run_newton_amap(centers, radii, rotations, c_ig, rho_ig, ...
                    N, newton_tol, max_newton_iter, cgm_tol, max_cgm_iter, verbose);

% Save iteration 1 data (for staged comparison)
h1 = history(1);
data = struct();
data.metadata = struct('case', 'annulus', 'N', N, 'm', m, ...
                       'stage', 'form_system', 'iteration', 1);
data.metadata.centers_re = real(centers).';
data.metadata.centers_im = imag(centers).';
data.metadata.radii = radii.';
data.metadata.rotations = rotations.';
data.matrices = struct('D', h1.D);
data.vectors = struct('g', h1.g);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'annulus_form_system_iter1.json'), data);
fprintf('  Wrote annulus_form_system_iter1.json\n');

data = struct();
data.metadata = struct('case', 'annulus', 'N', N, 'm', m, ...
                       'stage', 'solve_system', 'iteration', 1);
data.vectors = struct('U', h1.U);
data.scalars = struct();
data.matrices = struct();
save_reference_json(fullfile(output_dir, 'annulus_solve_system_iter1.json'), data);
fprintf('  Wrote annulus_solve_system_iter1.json\n');

data = struct();
data.metadata = struct('case', 'annulus', 'N', N, 'm', m, ...
                       'stage', 'newton_update', 'iteration', 1);
data.matrices = struct('S', h1.S, 'abs_eta', h1.abs_eta);
data.vectors = struct('c', h1.c, 'rho', h1.rho);
data.scalars = struct('normU', h1.normU);
save_reference_json(fullfile(output_dir, 'annulus_newton_update_iter1.json'), data);
fprintf('  Wrote annulus_newton_update_iter1.json\n');

% Save final converged state
data = struct();
data.metadata = struct('case', 'annulus', 'N', N, 'm', m, ...
                       'stage', 'converged', ...
                       'num_iterations', length(history));
data.matrices = struct('S', S_final, 'a', a_final);
data.vectors = struct('c', c_final, 'rho', rho_final);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'annulus_converged.json'), data);
fprintf('  Wrote annulus_converged.json (converged in %d iterations)\n', length(history));


%% ==================== Case 2: General m=3, N=64 ====================
fprintf('\n=== Generating general m=3 reference data ===\n');

N = 64;
centers = [0; 0.3+0.2i; -0.3-0.1i];
radii = [1.0; 0.1; 0.12];
rotations = [0; 0; 0];
c_ig = [0.3+0.2i; -0.3-0.1i];
rho_ig = [0.1; 0.12];
norm_cond = [1; 0; 0];
m = 3;

S0 = zeros(N, m);
for j = 1:m
    S0(:, j) = (0:N-1)' * 2*pi / N;
end

% --- P matrices ---
for nu = 1:m
    P = make_Pnu(nu, m, c_ig, rho_ig, N);
    data = struct();
    data.metadata = struct('case', 'general_m3', 'N', N, 'm', m, 'nu', nu, ...
                           'stage', 'P_matrix');
    data.matrices = struct('P', P);
    data.vectors = struct();
    data.scalars = struct();
    save_reference_json(fullfile(output_dir, ...
        sprintf('general_m3_P_nu%d.json', nu)), data);
    fprintf('  Wrote general_m3_P_nu%d.json\n', nu);
end

% --- Full Newton iteration ---
[S_final, c_final, rho_final, a_final, history] = ...
    run_newton_namap(centers, radii, rotations, c_ig, rho_ig, norm_cond, ...
                     N, newton_tol, max_newton_iter, cgm_tol, max_cgm_iter, verbose);

h1 = history(1);
data = struct();
data.metadata = struct('case', 'general_m3', 'N', N, 'm', m, ...
                       'stage', 'form_system', 'iteration', 1);
data.metadata.centers_re = real(centers).';
data.metadata.centers_im = imag(centers).';
data.metadata.radii = radii.';
data.metadata.rotations = rotations.';
data.metadata.norm_cond = norm_cond.';
data.matrices = struct('D', h1.D);
data.vectors = struct('g', h1.g);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'general_m3_form_system_iter1.json'), data);
fprintf('  Wrote general_m3_form_system_iter1.json\n');

data = struct();
data.metadata = struct('case', 'general_m3', 'N', N, 'm', m, ...
                       'stage', 'solve_system', 'iteration', 1);
data.vectors = struct('U', h1.U);
data.scalars = struct();
data.matrices = struct();
save_reference_json(fullfile(output_dir, 'general_m3_solve_system_iter1.json'), data);
fprintf('  Wrote general_m3_solve_system_iter1.json\n');

data = struct();
data.metadata = struct('case', 'general_m3', 'N', N, 'm', m, ...
                       'stage', 'newton_update', 'iteration', 1);
data.matrices = struct('S', h1.S, 'abs_eta', h1.abs_eta);
data.vectors = struct('c', h1.c, 'rho', h1.rho);
data.scalars = struct('normU', h1.normU);
save_reference_json(fullfile(output_dir, 'general_m3_newton_update_iter1.json'), data);
fprintf('  Wrote general_m3_newton_update_iter1.json\n');

data = struct();
data.metadata = struct('case', 'general_m3', 'N', N, 'm', m, ...
                       'stage', 'converged', ...
                       'num_iterations', length(history));
data.matrices = struct('S', S_final, 'a', a_final);
data.vectors = struct('c', c_final, 'rho', rho_final);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'general_m3_converged.json'), data);
fprintf('  Wrote general_m3_converged.json (converged in %d iterations)\n', length(history));


%% ==================== Case 3: Identity m=4 (th_gen_ex3), N=64 ====================
fprintf('\n=== Generating identity m=4 (th_gen_ex3) reference data ===\n');

N = 64;
centers = [0; -0.5; 0.25+0.43i; 0.25-0.43i];
radii = [1; 0.25; 0.25; 0.25];
rotations = [0; 0; 0; 0];
c_ig = [-0.4; 0.35+0.43i; 0.35-0.43i];
rho_ig = [0.25; 0.25; 0.25];
norm_cond = [1; 0; 0];
m = 4;

S0 = zeros(N, m);
for j = 1:m
    S0(:, j) = (0:N-1)' * 2*pi / N;
end

% --- P matrices ---
for nu = 1:m
    P = make_Pnu(nu, m, c_ig, rho_ig, N);
    data = struct();
    data.metadata = struct('case', 'identity_m4', 'N', N, 'm', m, 'nu', nu, ...
                           'stage', 'P_matrix');
    data.matrices = struct('P', P);
    data.vectors = struct();
    data.scalars = struct();
    save_reference_json(fullfile(output_dir, ...
        sprintf('identity_m4_P_nu%d.json', nu)), data);
    fprintf('  Wrote identity_m4_P_nu%d.json\n', nu);
end

% --- Full Newton iteration ---
[S_final, c_final, rho_final, a_final, history] = ...
    run_newton_namap(centers, radii, rotations, c_ig, rho_ig, norm_cond, ...
                     N, newton_tol, max_newton_iter, cgm_tol, max_cgm_iter, verbose);

h1 = history(1);
data = struct();
data.metadata = struct('case', 'identity_m4', 'N', N, 'm', m, ...
                       'stage', 'form_system', 'iteration', 1);
data.metadata.centers_re = real(centers).';
data.metadata.centers_im = imag(centers).';
data.metadata.radii = radii.';
data.metadata.rotations = rotations.';
data.metadata.norm_cond = norm_cond.';
data.matrices = struct('D', h1.D);
data.vectors = struct('g', h1.g);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'identity_m4_form_system_iter1.json'), data);
fprintf('  Wrote identity_m4_form_system_iter1.json\n');

data = struct();
data.metadata = struct('case', 'identity_m4', 'N', N, 'm', m, ...
                       'stage', 'solve_system', 'iteration', 1);
data.vectors = struct('U', h1.U);
data.scalars = struct();
data.matrices = struct();
save_reference_json(fullfile(output_dir, 'identity_m4_solve_system_iter1.json'), data);
fprintf('  Wrote identity_m4_solve_system_iter1.json\n');

data = struct();
data.metadata = struct('case', 'identity_m4', 'N', N, 'm', m, ...
                       'stage', 'newton_update', 'iteration', 1);
data.matrices = struct('S', h1.S, 'abs_eta', h1.abs_eta);
data.vectors = struct('c', h1.c, 'rho', h1.rho);
data.scalars = struct('normU', h1.normU);
save_reference_json(fullfile(output_dir, 'identity_m4_newton_update_iter1.json'), data);
fprintf('  Wrote identity_m4_newton_update_iter1.json\n');

data = struct();
data.metadata = struct('case', 'identity_m4', 'N', N, 'm', m, ...
                       'stage', 'converged', ...
                       'num_iterations', length(history));
data.matrices = struct('S', S_final, 'a', a_final);
data.vectors = struct('c', c_final, 'rho', rho_final);
data.scalars = struct();
save_reference_json(fullfile(output_dir, 'identity_m4_converged.json'), data);
fprintf('  Wrote identity_m4_converged.json (converged in %d iterations)\n', length(history));


fprintf('\n=== All reference data generated ===\n');
