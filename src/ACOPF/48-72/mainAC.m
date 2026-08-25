% ========================================================================
% Script: run_DistACOPF_all_modes.m   (UPDATED for generalized DistACOPF.m)
%
% Compatible with:
%   DistACOPF(partitionedDataFile, rho, residualThreshold, maxIterations, ...
%             centralizedCost, holds, mode, updateType, w1, w2, params)
%
% Runs DistACOPF for a specific partitioned case in:
%   1) Classic vs Primal-accelerated ADMM:
%        - Classic ('none')
%        - Primal Asymptotic    ('primal','asy')
%        - Primal Finite-time   ('primal','finite')
%        - Primal Fixed-time    ('primal','fixed')
%      -> Saves 2 plots: primal residuals & optimality gaps.
%
%   2) Classic vs Dual-accelerated ADMM:
%        - Classic ('none')
%        - Dual Finite-time     ('dual','finite')
%        - Dual Fixed-time      ('dual','fixed')
%      -> Saves 2 plots: primal residuals & optimality gaps.
%
% ========================================================================

clear; clc; close all;

%% -------------------------- User Parameters ----------------------------

% Example cases:
% pglib_opf_case48_ieee_2regions
% pglib_opf_case72_ieee_3regions
% pglib_opf_case118_ieee_3regions
% pglib_opf_case2869_pegase_12regions
caseName            = 'pglib_opf_case48_ieee_2regions';

partitionedDataFile = caseName;      % DistACOPF appends ".mat" if needed
loadedData = load(caseName);         % .mat file with "filename" inside
CasenameC  = loadedData.filename;    % MATPOWER .m case name

% ----------------------------------------------
%% CENTRALIZED AC-OPF (your helper name kept as-is)
% ----------------------------------------------
fprintf('> Running Centralized ACOPF...\n');
[resac, fvalac] = run_acopf_centralized(CasenameC); %#ok<ASGLU>

% ADMM/solver parameters
rho                = 1e0;
residualThreshold  = 1e-4;
maxIterations      = 100000;
holds              = 5e5;

% Centralized optimal cost (for optimality gap)
centralizedCost    = fvalac;     % set 0 to disable optimality gap plot

% Primal mode weights (weighted boundary copies)
w1_primal = 10.0;
w2_primal = 0.1;

% Output filenames (prefix only)
outPrefix = caseName;

%% ------------------- NEW: params struct for gains ----------------------
% You can tune these per your paper; defaults exist in DistACOPF, but we
% pass them explicitly to be consistent/reproducible.

params = struct();

% ---- Primal hat updates (asy/finite/fixed) ----
params.alpha1 = -1;     % asy
params.alpha2 = -0.1;     % finite/fixed
params.alpha3 = -0.05;     % fixed only
params.mu_p   = 0.50;      % 0 < mu < 1
params.eta_p  = 2.10;      % eta > 1

% ---- Dual updates (finite/fixed) ----
params.gamma1 = 1.0;     % finite/fixed
params.gamma2 = 1.0;     % fixed only
params.mu_d   = 0.10;      % 0 < mu < 1
params.eta_d  = -0.10;     % -1 < eta < 0
params.eps_d  = 1e-12;     % denom regularization

%% ======================= 1) PRIMAL MODES SET ===========================

fprintf('=== Running CLASSIC + PRIMAL modes for case %s ===\n', caseName);

%--- Classic ADMM (baseline) --------------------------------------------
% NOTE: updateType is ignored in mode='none', but we pass 'asy' for compatibility.
out_classic_primal = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                               maxIterations, centralizedCost, holds, ...
                               'none', 'asy', w1_primal, 0.0, params);

% --- Primal Asymptotic --------------------------------------------------
out_primal_asy = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                           maxIterations, centralizedCost, holds, ...
                           'primal', 'asy', w1_primal, w2_primal, params);

% --- Primal Finite-time -------------------------------------------------
out_primal_finite = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                              maxIterations, centralizedCost, holds, ...
                              'primal', 'finite', w1_primal, w2_primal, params);

% --- Primal Fixed-time --------------------------------------------------
out_primal_fixed = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                             maxIterations, centralizedCost, holds, ...
                             'primal', 'fixed', w1_primal, w2_primal, params);

%% --- Plot: Primal residual comparison (CLASSIC vs PRIMAL modes) --------
figure; hold on; grid on; box off;

k_c  = 1:out_classic_primal.iterationsRun;
k_pa = 1:out_primal_asy.iterationsRun;
k_pf = 1:out_primal_finite.iterationsRun;
k_px = 1:out_primal_fixed.iterationsRun;

semilogy(k_c,  out_classic_primal.errorLog,   'r:',  'LineWidth', 2.2);
semilogy(k_pa, out_primal_asy.errorLog,       'b--', 'LineWidth', 2.2);
semilogy(k_pf, out_primal_finite.errorLog,    'g-.', 'LineWidth', 2.2);
semilogy(k_px, out_primal_fixed.errorLog,     'k',   'LineWidth', 2.2);

% log ticks as 10^k
set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false));

xlabel('Iteration', 'FontSize', 20);
ylabel('Primal residual', 'FontSize', 20);
title(sprintf('%s', caseName), 'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg1 = legend({'Classic', 'Primal Asy', 'Primal Finite', 'Primal Fixed'}, 'Location','best');
set(leg1, 'FontSize', 16, 'Box', 'off');
set(gca, 'FontSize', 18, 'LineWidth', 1.5);

fname_primal_res = sprintf('%s_primal_residuals_ac.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_primal_res);

%% --- Plot: Optimality gap comparison (CLASSIC vs PRIMAL modes) ---------
figure; hold on; grid on; box off;

semilogy(k_c,  out_classic_primal.optimalityGap,   'r:',  'LineWidth', 2.2);
semilogy(k_pa, out_primal_asy.optimalityGap,       'b--', 'LineWidth', 2.2);
semilogy(k_pf, out_primal_finite.optimalityGap,    'g-.', 'LineWidth', 2.2);
semilogy(k_px, out_primal_fixed.optimalityGap,     'k',   'LineWidth', 2.2);

set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false));

xlabel('Iteration', 'FontSize', 20);
ylabel('Optimality gap [%]', 'FontSize', 20);
title(sprintf('%s', caseName), 'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg2 = legend({'Classic', 'Primal Asy', 'Primal Finite', 'Primal Fixed'}, 'Location','best');
set(leg2, 'FontSize', 16, 'Box', 'off');
set(gca, 'FontSize', 18, 'LineWidth', 1.5);

fname_primal_gap = sprintf('%s_primal_gap_ac.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_primal_gap);

fprintf('Saved PRIMAL comparison plots:\n  %s\n  %s\n', fname_primal_res, fname_primal_gap);

%% ======================== 2) DUAL MODES SET ============================

fprintf('=== Running CLASSIC + DUAL modes for case %s ===\n', caseName);

% --- Classic ADMM (baseline again) --------------------------------------
out_classic_dual = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                             maxIterations, centralizedCost, holds, ...
                             'none', 'asy', w1_primal, 0.0, params);

% --- Dual Finite-time ---------------------------------------------------
% NOTE: in 'dual' mode, w1/w2 are unused; pass 0 safely.
out_dual_finite = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                            maxIterations, centralizedCost, holds, ...
                            'dual', 'finite', 0.0, 0.0, params);

% --- Dual Fixed-time ----------------------------------------------------
out_dual_fixed = DistACOPF(partitionedDataFile, rho, residualThreshold, ...
                           maxIterations, centralizedCost, holds, ...
                           'dual', 'fixed', 0.0, 0.0, params);

%% --- Plot: Primal residual comparison (CLASSIC vs DUAL modes) ----------
figure; hold on; grid on; box off;

k_cd = 1:out_classic_dual.iterationsRun;
k_df = 1:out_dual_finite.iterationsRun;
k_dx = 1:out_dual_fixed.iterationsRun;

semilogy(k_cd, out_classic_dual.errorLog,   'r:',  'LineWidth', 2.2);
semilogy(k_df, out_dual_finite.errorLog,    'g-.', 'LineWidth', 2.2);
semilogy(k_dx, out_dual_fixed.errorLog,     'k',   'LineWidth', 2.2);

set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false));

xlabel('Iteration', 'FontSize', 20);
ylabel('Primal residual', 'FontSize', 20);
title(sprintf('%s', caseName), 'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg3 = legend({'Classic', 'Dual Finite', 'Dual Fixed'}, 'Location','best');
set(leg3, 'FontSize', 16, 'Box', 'off');
set(gca, 'FontSize', 18, 'LineWidth', 1.5);

fname_dual_res = sprintf('%s_dual_residuals_ac.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_dual_res);

%% --- Plot: Optimality gap comparison (CLASSIC vs DUAL modes) -----------
figure; hold on; grid on; box off;

semilogy(k_cd, out_classic_dual.optimalityGap,   'r:',  'LineWidth', 2.2);
semilogy(k_df, out_dual_finite.optimalityGap,    'g-.', 'LineWidth', 2.2);
semilogy(k_dx, out_dual_fixed.optimalityGap,     'k',   'LineWidth', 2.2);

set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, 'UniformOutput', false));

xlabel('Iteration', 'FontSize', 20);
ylabel('Optimality gap [%]', 'FontSize', 20);
title(sprintf('%s', caseName), 'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg4 = legend({'Classic', 'Dual Finite', 'Dual Fixed'}, 'Location','best');
set(leg4, 'FontSize', 16, 'Box', 'off');
set(gca, 'FontSize', 18, 'LineWidth', 1.5);

fname_dual_gap = sprintf('%s_dual_gap_ac.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_dual_gap);

fprintf('Saved DUAL comparison plots:\n  %s\n  %s\n', fname_dual_res, fname_dual_gap);

fprintf('=== All runs and plots completed for case %s ===\n', caseName);
