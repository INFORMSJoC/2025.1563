function out = DistDCOPF(partitionedDataFile, rho, residualThreshold, maxIterations, centralizedCost, holds, mode, updateType, w1, w2)
% =========================================================================
% Title       : Distributed DC Optimal Power Flow (DCOPF) Solver using ADMM
% Author      : Milad Hasanzadeh (extended with finite/fixed-time ADMM)
% Email       : e.mhasanzadeh1377@yahoo.com
% Affiliation : Department of Electrical and Computer Engineering,
%               Louisiana State University, Baton Rouge, LA, USA
% Date        : May 29, 2025  (extended version)
%
% Description :
% This function implements a distributed DC Optimal Power Flow (DCOPF) 
% solver based on the Alternating Direction Method of Multipliers (ADMM). 
% It now supports:
%
%   mode = 'none'   : classical ADMM (original behavior)
%   mode = 'primal' : primal-accelerated ADMM using auxiliary boundary
%                     angles \hat{\theta} with:
%                       updateType = 'asy'    (asymptotic)
%                       updateType = 'finite' (finite-time)
%                       updateType = 'fixed'  (fixed-time)
%   mode = 'dual'   : dual-accelerated ADMM using finite/fixed-time
%                     update rules for the dual variables:
%                       updateType = 'finite'
%                       updateType = 'fixed'
%
% Optional:
%   w1, w2 : weighting coefficients in the weighted consistency residual
%            for the copied boundary angles, used only in mode='primal'.
%            If omitted, defaults are w1 = 1, w2 = 0 (classical ADMM).
%
% FUNCTION SIGNATURE:
%   out = DistDCOPF(partitionedDataFile, rho, residualThreshold, ...
%                   maxIterations, centralizedCost, holds, ...
%                   mode, updateType, w1, w2)
%
% INPUTS:
%   partitionedDataFile : string, name of the partitioned .mat file
%   rho                 : ADMM penalty parameter (for augmented Lagrangian)
%   residualThreshold   : stopping threshold on worst primal residual
%   maxIterations       : maximum number of ADMM iterations
%   centralizedCost     : centralized optimal total cost (dollars); 
%                         set to 0 to skip optimality gap computation
%   holds               : upper bound on how far rho can adapt
%
%   mode                : 'none' | 'primal' | 'dual'
%   updateType          : 
%                         mode='none'  -> ignored (classical ADMM)
%                         mode='primal'-> 'asy' | 'finite' | 'fixed'
%                         mode='dual'  -> 'finite' | 'fixed'
%   w1, w2              : optional weights for primal mode
%
% OUTPUT:
%   out : struct containing key results:
%       .converged          - logical, true if converged before maxIterations
%       .iterations         - number of iterations actually run
%       .errorLog           - vector of max error per iteration
%       .optimalityGap      - vector of optimality gap per iteration (NaN if centralizedCost==0)
%       .lambda             - final dual variables (2*nit x 1)
%       .currentObjective   - objective value at final iteration
%       .tieLines           - tie-line region pairing cell array
%       .total_Pd           - total demand (MW)
%       .total_Pg           - total generation (MW)
%       .nit                - total number of interfaces
%       .numRegions         - number of regions
% =========================================================================

%% -------------------------- Handle Inputs ------------------------------
if ~endsWith(partitionedDataFile, '.mat')
    partitionedDataFile = [partitionedDataFile, '.mat'];
end

if ~isfile(partitionedDataFile)
    error('Required file "%s" not found. Please ensure it is in the current directory.', ...
        partitionedDataFile);
end

if nargin < 5 || isempty(centralizedCost)
    centralizedCost = 0;
end
if nargin < 6 || isempty(holds)
    holds = 1e8;
end

% Mode and inner update type
if nargin < 7 || isempty(mode)
    mode = 'none';       % classical ADMM
end
if nargin < 8 || isempty(updateType)
    updateType = 'asy';  % default if mode='primal' and not specified
end
mode       = lower(mode);
updateType = lower(updateType);

% Optional weights for primal acceleration
if nargin < 9 || isempty(w1)
    w1 = 1;              % classical ADMM recovered if w2 = 0
end
if nargin < 10 || isempty(w2)
    w2 = 0;              % no auxiliary contribution by default
end

% Flags
usePrimalAccel = strcmpi(mode, 'primal');
useDualAccel   = strcmpi(mode, 'dual');

%% ---------------------- Acceleration Parameters ------------------------
% These should satisfy the theoretical ranges from the paper.
% You can tune them for performance.



% --- Primal auxiliary update gains (for \hat{\theta}) ---
alpha1    = -1;   % asymptotic gain (alpha_1 < 0)
alpha2    = -0.1;   % finite/fixed gain (alpha_2 < 0)
alpha3    = -0.05;   % fixed-time extra gain (alpha_3 < 0)
muPrimal  = 0.5;    % 0 < mu < 1
etaPrimal = 2.1;    % > 1

% --- Dual gradient-shaping gains (for lambda) ---

muDual    = 0.1;    % 0 < mu < 1
etaDual   = -0.1;   % -1 < eta < 0

gamma1    = .5;   % finite/fixed gain (gamma_1 < 0)
gamma2    = 1.2;   % fixed-time extra gain (gamma_2 < 0)

%% ----------------------------- Load Case -------------------------------
loadedData   = load(partitionedDataFile);
numRegions   = loadedData.num_regions;
totalInterfaces = 0;

%% ------------- Augment Each Region with Interregional Tie-Lines --------
for regionIdx = 1:numRegions
    regionID      = ['R' num2str(regionIdx)];
    regionName    = ['mpc_region', regionID];
    tielinesName  = ['interregional_tielines', regionID];
    baseBranch    = loadedData.(regionName).branch;

    % Convert global interregional table to cell
    dataMatrix = table2cell(loadedData.interregional_tielines_total);

    % Rows where this region participates in the tie-lines
    row_indices = strcmp(dataMatrix(:, 1), regionName) | ...
                  strcmp(dataMatrix(:, 2), regionName);

    filteredData = dataMatrix(row_indices, :);
    nbus        = size(loadedData.(regionName).bus, 1);
    nTieLocal   = size(filteredData, 1);

    % Safety: handle regions with no tie-lines
    if nTieLocal == 0
        loadedData.(tielinesName) = cell(0, size(filteredData, 2));
        continue;
    end

    % Initialize check matrix per region (direction indicator)
    checkmattrix = zeros(nTieLocal, 1);

    for i = 1:nTieLocal
        region1_str = filteredData{i, 1};
        region2_str = filteredData{i, 2};

        region1_num = sscanf(region1_str, 'mpc_regionR%d');
        region2_num = sscanf(region2_str, 'mpc_regionR%d');

        if region1_num == regionIdx
            checkmattrix(i) = 1;
        else
            checkmattrix(i) = -1;
        end
    end

    % Numeric tie-line rows (branch columns)
    output_rows = cell2mat(filteredData(:, 3:end));

    % Adjust from/to bus indices for out-of-region buses
    for i = 1:nTieLocal
        if checkmattrix(i) == 1
            % This region is "region1": to-bus is external
            output_rows(i, 2) = nbus + i;
        else
            % This region is "region2": from-bus is external
            output_rows(i, 1) = nbus + i;
        end
    end

    % Append tie-lines to region branch matrix
    tieLines = output_rows;
    loadedData.(regionName).branch = [baseBranch; tieLines];

    % Store the numeric tie-line row in column 3 for this region’s tie-line cell
    for i = 1:nTieLocal
        filteredData{i, 3} = output_rows(i, :);  % 1×13 double
    end

    % Save region-specific tie-line data
    loadedData.(tielinesName) = filteredData;
end

%% ---------------------- Initialize Region Variables --------------------
for regionIdx = 1:numRegions
    regionID      = ['R' num2str(regionIdx)];
    regionName    = ['mpc_region', regionID];
    tielinesName  = ['interregional_tielines', regionID];
    regionData    = loadedData.(regionName);
    tielinesData  = loadedData.(tielinesName);

    % Create variable names
    nbName        = ['nb', regionID];
    nitName       = ['nit', regionID];
    conName       = ['con', regionID];
    yonName       = ['Yon', regionID];
    phaseName     = ['phase_angles_', regionID];
    hatPhaseName  = ['hat_phase_angles_', regionID];   % NEW: auxiliary angles
    regionVar     = ['region_', regionID];
    RealpowerName = ['Realpowerflow_', regionID];

    % Extract region data
    nb       = size(regionData.bus, 1);
    nitLocal = size(tielinesData, 1);
    branchMatrix = regionData.branch;

    % Initialize values in base workspace (for DCOPF_SP script)
    assignin('base', nbName, nb);
    assignin('base', nitName, nitLocal);
    assignin('base', phaseName, zeros(2 * nitLocal, 1));
    assignin('base', RealpowerName, zeros(nitLocal, 1));
    assignin('base', conName, []);
    assignin('base', yonName, []);

    % NEW: initialize auxiliary boundary angles for primal mode
    assignin('base', hatPhaseName, ones(nitLocal, 1));

    % Extract tie-line rows (buses connected to this region)
    outRegionRows = branchMatrix(:, 1) > nb | branchMatrix(:, 2) > nb;
    tieLineData   = branchMatrix(outRegionRows, :);

    regionBusIDs  = cell(size(tieLineData, 1), 1);
    for j = 1:size(tieLineData, 1)
        if tieLineData(j, 1) <= nb
            regionBusIDs{j} = tieLineData(j, 1);
        elseif tieLineData(j, 2) <= nb
            regionBusIDs{j} = tieLineData(j, 2);
        end
    end
    regionBusIDs = regionBusIDs(~cellfun('isempty', regionBusIDs));
    assignin('base', regionVar, regionBusIDs);

    totalInterfaces = totalInterfaces + nitLocal;
end

% Total number of interfaces (each interface appears twice across regions)
nit = totalInterfaces / 2;
assignin('base', 'nit', nit);

% Load original MATPOWER case only to get baseMVA
mpc0   = loadcase(loadedData.filename);
nb_full = size(mpc0.bus, 1);
baseMW  = mpc0.baseMVA;

% ---- Angle scaling for normalized residuals ----
x      = mpc0.branch(:, 4);
rateA  = mpc0.branch(:, 6);
theta_max_approx = abs(x) .* (rateA / baseMW);

valid = theta_max_approx > 0 & ~isnan(theta_max_approx) & ~isinf(theta_max_approx);
if any(valid)
    theta_scale = median(theta_max_approx(valid));   % e.g., ~0.1–0.3 rad typically
else
    theta_scale = 0.1;   % safe fallback if no valid limits
end
theta_scale = 1;
%% ------------------------ Get ADMM Parameters -------------------------
lambda        = zeros(2 * nit, 1);
errorValues   = zeros(2 * nit, 1);
errorLog      = zeros(maxIterations, 1);
optimalityGap = NaN(maxIterations, 1);
grad          = zeros(nit, 1);
gradDual      = zeros(nit, 1);
currentObjective = 0;

% Residual vectors for adaptive rho
r_theta      = zeros(2 * nit, 1);   % current normalized residuals
r_theta_old  = zeros(2 * nit, 1);   % previous iteration residuals

% Heuristic parameters for adaptive rho
mu        = 10;   % residual balance factor
tau_incr  = 1.5;    % how much to increase rho
tau_decr  = 1.5;   % how much to decrease rho
eta       = 1;    % ADMM dual step-size scaling (not the exponent eta in the paper!)
max_count = 1;


rhoLog  = zeros(maxIterations, 1);   % to store rho per iteration
dualLog = zeros(maxIterations, 1);
iter_count = 0;
Flag_rho   = 0;

%% ---------------- Prepare Interregional Interface Buses ----------------
tieLineTable = loadedData.interregional_tielines_total;
tieLines     = table2cell(tieLineTable);

% Reorder region strings so that region1_num < region2_num
for i = 1:nit
    region1_str = tieLines{i, 1};
    region2_str = tieLines{i, 2};

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    if region1_num > region2_num
        [tieLines{i, 1}, tieLines{i, 2}] = deal(region2_str, region1_str);
    end
end

% Replace with sorted ordering
tieLineTable = cell2table(tieLines);

% Precompute interconnection counters (local interface index per region)
interfaceCounter = zeros(numRegions, 1);

% NEW: interface map (global -> (region1, localIdx1; region2, localIdx2))
interfaceMap = struct('r1', cell(nit, 1), ...
                      'r2', cell(nit, 1), ...
                      'idx1', cell(nit, 1), ...
                      'idx2', cell(nit, 1));

for i = 1:nit
    region1_str = tieLines{i, 1};  % e.g., 'mpc_regionR3'
    region2_str = tieLines{i, 2};  % e.g., 'mpc_regionR5'

    region1_num = sscanf(region1_str, 'mpc_regionR%d');
    region2_num = sscanf(region2_str, 'mpc_regionR%d');

    for j = 1:numRegions
        if j == region1_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            a = interfaceCounter(j);
        elseif j == region2_num
            interfaceCounter(j) = interfaceCounter(j) + 1;
            b = interfaceCounter(j);
        end
    end

    % Store mapping between global interface i and local indices in two regions
    interfaceMap(i).r1   = region1_num;
    interfaceMap(i).r2   = region2_num;
    interfaceMap(i).idx1 = a;
    interfaceMap(i).idx2 = b;

    evalin('base', sprintf('conR%d = [conR%d; %d];', region1_num, region1_num, b));
    evalin('base', sprintf('conR%d = [conR%d; %d];', region2_num, region2_num, a));
    evalin('base', sprintf('YonR%d = [YonR%d; %d];', region1_num, region1_num, i));
    evalin('base', sprintf('YonR%d = [YonR%d; %d];', region2_num, region2_num, i));
end
totalSubproblemWallTime = 0;
%% -------------------------- Begin ADMM Loop ----------------------------
fprintf('Starting ADMM iterations...\n');

for iter = 1:maxIterations

    %% --- NEW: Primal auxiliary update (if mode = 'primal') -------------
    if usePrimalAccel
        % Update \hat{\theta}_j^{(r)} for all interfaces based on the
        % previously stored boundary angles (phase_angles_R*).
        for iif = 1:nit
            r1   = interfaceMap(iif).r1;
            r2   = interfaceMap(iif).r2;
            idx1 = interfaceMap(iif).idx1;
            idx2 = interfaceMap(iif).idx2;

            reg1 = ['R' num2str(r1)];
            reg2 = ['R' num2str(r2)];

            % Owned bus angles on each side (from previous iteration)
            angle1 = evalin('base', ['phase_angles_' reg1]);
            angle2 = evalin('base', ['phase_angles_' reg2]);

            if isempty(angle1) || isempty(angle2)
                % Should not happen, but just in case
                continue;
            end

            theta1_owned = angle1(idx1);
            theta2_owned = angle2(idx2);

            hat1 = evalin('base', ['hat_phase_angles_' reg1]);
            hat2 = evalin('base', ['hat_phase_angles_' reg2]);

            % Asymptotic / finite / fixed-time updates
            if strcmpi(updateType, 'asy') || strcmpi(updateType, 'asym')
                e1 = hat1(idx1) - theta2_owned;
                e2 = hat2(idx2) - theta1_owned;

                hat1(idx1) = hat1(idx1) + alpha1 * e1;
                hat2(idx2) = hat2(idx2) + alpha1 * e2;

            elseif strcmpi(updateType, 'finite')
                e1 = hat1(idx1) - theta2_owned;
                e2 = hat2(idx2) - theta1_owned;

                hat1(idx1) = hat1(idx1) + alpha2 * sign(e1) * abs(e1)^muPrimal;
                hat2(idx2) = hat2(idx2) + alpha2 * sign(e2) * abs(e2)^muPrimal;

            elseif strcmpi(updateType, 'fixed')
                e1 = hat1(idx1) - theta2_owned;
                e2 = hat2(idx2) - theta1_owned;

                hat1(idx1) = hat1(idx1) + ...
                    alpha3 * sign(e1) * abs(e1)^muPrimal + ...
                    alpha3 * sign(e1) * abs(e1)^etaPrimal;

                hat2(idx2) = hat2(idx2) + ...
                    alpha3 * sign(e2) * abs(e2)^muPrimal + ...
                    alpha3 * sign(e2) * abs(e2)^etaPrimal;
            end

            assignin('base', ['hat_phase_angles_' reg1], hat1);
            assignin('base', ['hat_phase_angles_' reg2], hat2);
        end
    end

    %% --- Local DCOPF solve for each region ----------------------------
    currentObjective = 0;
    Pg_total = 0;

    for r = 1:numRegions
        reg = ['R' num2str(r)];

        % DCOPF_SP is a script that uses:
        %  - loadedData, reg, lambda, rho, nit, theta_scale, mode, updateType, w1, w2
        %  - plus region-specific variables from 'base' (nbR#, nitR#, conR#, YonR#, region_R#)
        run('DCOPF_SP.m');
       
        totalSubproblemWallTime = totalSubproblemWallTime + local_solver_time;


        nitVal        = evalin('base', ['nit', reg]);
        nbVal         = evalin('base', ['nb', reg]);
        regionBus     = evalin('base', ['region_', reg]);
        phaseName     = ['phase_angles_', reg];
        RealpowerName = ['Realpowerflow_', reg];

        % Reinitialize per region to avoid leftover values
        angle  = zeros(2 * nitVal, 1);
        Real_p = zeros(nitVal, 1);

        % Extract updated angles and "flow" proxy
        for k = 1:nitVal
            angle(k, 1)          = value(Delta(regionBus{k}));
            angle(k + nitVal, 1) = value(Delta(k + nbVal));
            Real_p(k, 1)         = value(Delta(regionBus{k})) - ...
                                   value(Delta(k + nbVal));
        end
        assignin('base', phaseName, angle);
        assignin('base', RealpowerName, Real_p);

        % Local objective using MATPOWER's totcost
        Pg_MW  = value(Pg) * baseMW;      % column vector [ngen_r x 1]
        cost_vec = totcost(mpc.gencost, Pg_MW);
        cost_reg = sum(cost_vec(:));

        currentObjective = currentObjective + cost_reg;
        Pg_total         = Pg_total + sum(value(Pg));
    end

    %% --- Dual Update Step ----------------------------------------------
    interfaceCounter = zeros(numRegions, 1);  % reset per iteration

    for i = 1:nit
        str1 = tieLines{i, 1};  % 'mpc_regionR#'
        str2 = tieLines{i, 2};

        region1_num = sscanf(str1, 'mpc_regionR%d');
        region2_num = sscanf(str2, 'mpc_regionR%d');

        for j = 1:numRegions
            if j == region1_num
                interfaceCounter(j) = interfaceCounter(j) + 1;
                a = interfaceCounter(j);
            elseif j == region2_num
                interfaceCounter(j) = interfaceCounter(j) + 1;
                b = interfaceCounter(j);
            end
        end

        reg1 = ['R' num2str(region1_num)];
        reg2 = ['R' num2str(region2_num)];

    angle1 = evalin('base', ['phase_angles_' reg1]);
    angle2 = evalin('base', ['phase_angles_' reg2]);
    nit1   = evalin('base', ['nit' reg1]);
    nit2   = evalin('base', ['nit' reg2]);

    % ------------------ Residuals depend on mode -----------------------
    % if usePrimalAccel
    % 
    % 
    %     hat1 = evalin('base', ['hat_phase_angles_' reg1]);  % length nit1
    %     hat2 = evalin('base', ['hat_phase_angles_' reg2]);  % length nit2
    % 
    %     denom = w1 + w2;
    %     if denom == 0
    %         denom = 1;   % fallback: just use tilde_theta
    %     end
    % 
    %     eff_copy_i_from_reg2 = (w1 * angle2(nit2 + b) + w2 * hat2(b)) / denom;
    %     grad(i)              = angle1(a) - eff_copy_i_from_reg2;
    % 
    %     eff_copy_j_from_reg1 = (w1 * angle1(nit1 + a) + w2 * hat1(a)) / denom;
    %     gradDual(i)          = eff_copy_j_from_reg1 - angle2(b);
    % 
    % else
        % ---------------- Classic / dual modes -------------------------
    
        grad(i)     = angle1(a)        - angle2(nit2 + b);
        gradDual(i) = angle1(nit1 + a) - angle2(b);
    % end

    % ---------------- Normalize residuals -------------------------------
    r_norm      = grad(i)     / theta_scale;
    rdual_norm  = gradDual(i) / theta_scale;

    % Store normalized residuals in r_theta (length 2*nit)
    r_theta(i)      = r_norm;
    r_theta(i+nit)  = rdual_norm;

    % Track absolute normalized errors
    errorValues(i)      = abs(r_norm);
    errorValues(i+nit)  = abs(rdual_norm);

    % ---------------- Lambda update (unchanged logic) -------------------
    if ~useDualAccel
        % Classical ADMM dual update
        lambda(i)      = lambda(i)      + eta * rho * r_norm;
        lambda(i+nit)  = lambda(i+nit)  + eta * rho * rdual_norm;
    else
        % Dual-accelerated finite / fixed-time update on the
        % normalized gradient r_theta.
        g1 = r_norm;
        g2 = rdual_norm;rhobar = rho/(rho^muDual);

        if strcmpi(updateType, 'finite')
            % Finite-time
            if abs(g1) > 0
                lambda(i) = lambda(i) + ...
                    rhobar *gamma1 * sign(g1) * abs(g1)^(1 - muDual);
            end
            if abs(g2) > 0
                lambda(i+nit) = lambda(i+nit) + ...
                    rhobar *gamma1 * sign(g2) * abs(g2)^(1 - muDual);
            end

        elseif strcmpi(updateType, 'fixed')
            % Fixed-time
            if abs(g1) > 0
                lambda(i) = lambda(i) + ...
                    rhobar *gamma2 * sign(g1) * abs(g1)^(1 - muDual) + ...
                    rhobar *gamma2 * sign(g1) * abs(g1)^(1 - etaDual);
            end
            if abs(g2) > 0
                lambda(i+nit) = lambda(i+nit) + ...
                    rhobar *gamma2 * sign(g2) * abs(g2)^(1 - muDual) + ...
                    rhobar *gamma2 * sign(g2) * abs(g2)^(1 - etaDual);
            end
        else
            % Fallback: classical ADMM if invalid updateType
            lambda(i)      = lambda(i)      + eta * rho * r_norm;
            lambda(i+nit)  = lambda(i+nit)  + eta * rho * rdual_norm;
        end
    end

    end

    %% --- Log and Adaptive rho ------------------------------------------
    % Primal residual = max normalized consensus mismatch
    primalRes       = max(abs(r_theta));
    errorLog(iter)  = primalRes;

    % Dual residual = change in residuals between iterations
    dualRes         = norm(r_theta - r_theta_old, inf);
    dualLog(iter)   = dualRes;
    rhoLog(iter)    = rho;

    if Flag_rho == 0
        % Boyd-style adaptive rho (still uses primal/dual residuals)
        if primalRes > mu * dualRes
            if rho < holds
            rho     = rho * tau_incr; Flag_rho = 1;
            end
        elseif dualRes > mu * primalRes
            rho     = rho / tau_decr; Flag_rho = 1;
        end
    else
        iter_count = iter_count + 1;
    end

    if iter_count == max_count
        iter_count = 1; Flag_rho = 0;
    end

    % Store residual for next iteration
    r_theta_old = r_theta;

    if centralizedCost ~= 0
        optimalityGap(iter) = abs(((currentObjective) - centralizedCost) ...
                                   / centralizedCost) * 100;
        fprintf('Iter %3d | Error: %.6e | Obj: %.4f | Gap: %.3f %% | rho: %.4f | mode: %s/%s\n', ...
            iter, primalRes, currentObjective, optimalityGap(iter), rho, mode, updateType);
    else
        fprintf('Iter %3d | Error: %.6e | Obj: %.4f | rho: %.4f | mode: %s/%s\n', ...
            iter, primalRes, currentObjective, rho, mode, updateType);
    end

    % Convergence check
    if primalRes < residualThreshold
        fprintf('Convergence achieved at iteration %d\n', iter);
        break;
    end
end

converged = (iter < maxIterations);

if ~converged
    disp('Maximum iterations reached without full convergence.');
end

% Report total wall time spent solving all local DCOPF subproblems
fprintf('Total sum of DCOPF subproblem solver time: %.4f seconds\n', ...
    totalSubproblemWallTime);


%% ---------------------- Post-Processing: Demand/Generation ------------
mpc = loadcase(loadedData.filename);

Pd = mpc.bus(:, 3);   % Real power demand (MW)
total_Pd = sum(Pd);   % Total real power demand (MW)

total_Pg = (mpc.baseMVA) .* Pg_total;

fprintf('Total Real Power Demand (Pd): %.2f MW\n', total_Pd);
fprintf('Total Real Power Generation (Pg): %.2f MW\n', total_Pg);

%% ---------------------------- Build Output Struct ----------------------
out = struct();
out.converged        = converged;
out.iterations       = iter;
out.errorLog         = errorLog(1:iter);
out.optimalityGap    = optimalityGap(1:iter);
out.lambda           = lambda;
out.currentObjective = currentObjective;
out.tieLines         = tieLines;
out.total_Pd         = total_Pd;
out.total_Pg         = total_Pg;
out.nit              = nit;
out.numRegions       = numRegions;
% NEW: total subproblem wall time over the whole ADMM run
out.totalSubproblemWallTime = totalSubproblemWallTime;
end
