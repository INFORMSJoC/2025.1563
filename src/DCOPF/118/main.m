% ========================================================================
% Script: run_DistDCOPF_all_modes.m
%
% Runs DistDCOPF for a specific partitioned case in:
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
% Assumes the extended DistDCOPF.m and DCOPF_SP.m with mode/updateType.
% ========================================================================

clear; clc; close all;

%% -------------------------- User Parameters ----------------------------

% pglib_opf_case48_ieee_2regions pglib_opf_case72_ieee_3regions
% pglib_opf_case118_ieee_3regions pglib_opf_case2869_pegase_12regions

caseName            = 'pglib_opf_case118_ieee_3regions';  

% If needed, automatic .mat extension happens inside DistDCOPF.
partitionedDataFile = caseName;
loadedData = load(caseName);        % .mat file with "filename" inside
CasenameC  = loadedData.filename;   % MATPOWER .m case name

% ----------------------------------------------
%% CENTRALIZED DC-OPF
%% ----------------------------------------------
fprintf('> Running Centralized DCOPF...\n');
[resdc, fvaldc] = run_dcopf_centralized(CasenameC);

% ADMM/solver parameters
rho                = 1e2;      % base ADMM penalty (same as your main setup)
residualThreshold  = 1e-4;     % stopping criterion on max normalized residual
maxIterations      = 100000;   % cap on ADMM iterations
holds              = 2e5;      % upper bound for adaptive rho

% >>> Centralized optimal cost (for optimality gap).
%     Set this to your known centralized DC OPF cost for this case.
%     If you set centralizedCost = 0, the optimality gap curves will be NaN.
centralizedCost    =  fvaldc;   % <-- EDIT THIS to your centralized cost

% Primal mode weights (for weighted boundary copies)
w1_primal = 10.0;
w2_primal = 0.1;   % you can tune this; w2=0 recovers classical copy

% Output filenames (prefix only; extension added automatically)
outPrefix = caseName;   % you can add more formatting if you want

%% ======================= 1) PRIMAL MODES SET ===========================

fprintf('=== Running CLASSIC + PRIMAL modes for case %s ===\n', caseName);

% --- Classic ADMM (baseline) --------------------------------------------
out_classic_primal = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                               maxIterations, centralizedCost, holds, ...
                               'none', 'asy', w1_primal, 0.0);

% --- Primal Asymptotic --------------------------------------------------
out_primal_asy = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                           maxIterations, centralizedCost, holds, ...
                           'primal', 'asy', w1_primal, w2_primal);

% --- Primal Finite-time -------------------------------------------------
out_primal_finite = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                              maxIterations, centralizedCost, holds, ...
                              'primal', 'finite', w1_primal, w2_primal);

% --- Primal Fixed-time --------------------------------------------------
out_primal_fixed = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                             maxIterations, centralizedCost, holds, ...
                             'primal', 'fixed', w1_primal, w2_primal);

%% --- Plot: Primal residual comparison (CLASSIC vs PRIMAL modes) --------
figure;
hold on; grid on; box off;

k_c  = 1:out_classic_primal.iterations;
k_pa = 1:out_primal_asy.iterations;
k_pf = 1:out_primal_finite.iterations;
k_px = 1:out_primal_fixed.iterations;

% Thicker lines for publication
semilogy(k_c,  out_classic_primal.errorLog,   'r:',  'LineWidth', 2.2);  % classic
semilogy(k_pa, out_primal_asy.errorLog,       'b--', 'LineWidth', 2.2);  % primal asy
semilogy(k_pf, out_primal_finite.errorLog,    'g-.', 'LineWidth', 2.2);  % primal finite
semilogy(k_px, out_primal_fixed.errorLog,     'k',   'LineWidth', 2.2);  % primal fixed

% --- Log-scale ticks as 10^k -------------------------------------------
set(gca, 'YScale', 'log');  % already implied by semilogy but explicit
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, ...
                     'UniformOutput', false));
% -----------------------------------------------------------------------

% Larger font sizes for publication
xlabel('Iteration', ...
    'FontSize', 20, 'Interpreter', 'none');
ylabel('Primal residual', ...
    'FontSize', 20, 'Interpreter', 'none');
title(sprintf('%s', caseName), ...
      'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg1 = legend({'Classic', 'Primal Asy', 'Primal Finite', 'Primal Fixed'}, ...
              'Location','best');
set(leg1, 'FontSize', 16, 'Box', 'off'); 

set(gca, 'FontSize', 18, ...          % axis tick labels
         'LineWidth', 1.5);           % axis box line width

% Save figure
fname_primal_res = sprintf('%s_primal_residuals_dc.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');   % better export
saveas(gcf, fname_primal_res);

%% --- Plot: Optimality gap comparison (CLASSIC vs PRIMAL modes) ---------
figure;
hold on; grid on; box off;

semilogy(k_c,  out_classic_primal.optimalityGap,   'r:',  'LineWidth', 2.2);
semilogy(k_pa, out_primal_asy.optimalityGap,       'b--', 'LineWidth', 2.2);
semilogy(k_pf, out_primal_finite.optimalityGap,    'g-.', 'LineWidth', 2.2);
semilogy(k_px, out_primal_fixed.optimalityGap,     'k',   'LineWidth', 2.2);

% --- Log-scale ticks as 10^k -------------------------------------------
set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, ...
                     'UniformOutput', false));
% -----------------------------------------------------------------------

xlabel('Iteration', ...
    'FontSize', 20, 'Interpreter', 'none');
ylabel('Optimality gap [%]', ...
    'FontSize', 20, 'Interpreter', 'none');
title(sprintf('%s', caseName), ...
      'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg2 = legend({'Classic', 'Primal Asy', 'Primal Finite', 'Primal Fixed'}, ...
              'Location','best');
set(leg2, 'FontSize', 16, 'Box', 'off');

set(gca, 'FontSize', 18, ...
         'LineWidth', 1.5);

% Save figure
fname_primal_gap = sprintf('%s_primal_gap_dc.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_primal_gap);

fprintf('Saved PRIMAL comparison plots:\n  %s\n  %s\n', ...
        fname_primal_res, fname_primal_gap);

%% ======================== 2) DUAL MODES SET ============================

fprintf('=== Running CLASSIC + DUAL modes for case %s ===\n', caseName);

% --- Classic ADMM (baseline again, for dual comparison) -----------------
out_classic_dual = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                             maxIterations, centralizedCost, holds, ...
                             'none', 'asy', w1_primal, 0.0);

% --- Dual Finite-time ---------------------------------------------------
out_dual_finite = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                            maxIterations, centralizedCost, holds, ...
                            'dual', 'finite', w1_primal, 0.0);

% --- Dual Fixed-time ----------------------------------------------------
out_dual_fixed = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
                           maxIterations, centralizedCost, holds, ...
                           'dual', 'fixed', w1_primal, 0.0);

%% --- Plot: Primal residual comparison (CLASSIC vs DUAL modes) ----------
figure;
hold on; grid on; box off;

k_cd = 1:out_classic_dual.iterations;
k_df = 1:out_dual_finite.iterations;
k_dx = 1:out_dual_fixed.iterations;

semilogy(k_cd, out_classic_dual.errorLog,   'r:',  'LineWidth', 2.2);  % classic
semilogy(k_df, out_dual_finite.errorLog,    'g-.', 'LineWidth', 2.2);  % dual finite
semilogy(k_dx, out_dual_fixed.errorLog,     'k',   'LineWidth', 2.2);  % dual fixed

% --- Log-scale ticks as 10^k -------------------------------------------
set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, ...
                     'UniformOutput', false));
% -----------------------------------------------------------------------

xlabel('Iteration', ...
    'FontSize', 20, 'Interpreter', 'none');
ylabel('Primal residual', ...
    'FontSize', 20, 'Interpreter', 'none');
title(sprintf('%s', caseName), ...
      'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg3 = legend({'Classic', 'Dual Finite', 'Dual Fixed'}, 'Location','best');
set(leg3, 'FontSize', 16, 'Box', 'off');

set(gca, 'FontSize', 18, ...
         'LineWidth', 1.5);

% Save figure
fname_dual_res = sprintf('%s_dual_residuals_dc.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_dual_res);

%% --- Plot: Optimality gap comparison (CLASSIC vs DUAL modes) -----------
figure;
hold on; grid on; box off;

semilogy(k_cd, out_classic_dual.optimalityGap,   'r:',  'LineWidth', 2.2);
semilogy(k_df, out_dual_finite.optimalityGap,    'g-.', 'LineWidth', 2.2);
semilogy(k_dx, out_dual_fixed.optimalityGap,     'k',   'LineWidth', 2.2);

% --- Log-scale ticks as 10^k -------------------------------------------
set(gca, 'YScale', 'log');
yl = ylim;
exp_min = floor(log10(yl(1)));
exp_max = ceil(log10(yl(2)));
exponents = exp_min:exp_max;
yticks(10.^exponents);
yticklabels(arrayfun(@(e) sprintf('10^{%d}', e), exponents, ...
                     'UniformOutput', false));
% -----------------------------------------------------------------------

xlabel('Iteration', ...
    'FontSize', 20, 'Interpreter', 'none');
ylabel('Optimality gap [%]', ...
    'FontSize', 20, 'Interpreter', 'none');
title(sprintf('%s', caseName), ...
      'Interpreter','none', 'FontSize', 22, 'FontWeight', 'bold');

leg4 = legend({'Classic', 'Dual Finite', 'Dual Fixed'}, 'Location','best');
set(leg4, 'FontSize', 16, 'Box', 'off');

set(gca, 'FontSize', 18, ...
         'LineWidth', 1.5);

% Save figure
fname_dual_gap = sprintf('%s_dual_gap_dc.png', outPrefix);
set(gcf, 'PaperPositionMode', 'auto');
saveas(gcf, fname_dual_gap);

fprintf('Saved DUAL comparison plots:\n  %s\n  %s\n', ...
        fname_dual_res, fname_dual_gap);

fprintf('=== All runs and plots completed for case %s ===\n', caseName);
