function out = DistACOPF(partitionedDataFile, rho, residualThreshold, maxIterations, centralizedCost, holds, mode, updateType, w1, w2, params)
% =========================================================================
% Generalized Distributed ACOPF (ADMM) Solver
% Supports:
%   mode = 'none'   : classical ADMM (no auxiliary hats, classical dual updates)
%   mode = 'primal' : primal-accelerated ADMM using auxiliary boundary voltages (hat)
%                    updateType = 'asy' | 'finite' | 'fixed'
%                    weighted consistency uses w1,w2
%   mode = 'dual'   : dual-accelerated ADMM (finite/fixed-time dual updates)
%                    updateType = 'finite' | 'fixed'
%
% Inputs:
%   partitionedDataFile, rho, residualThreshold, maxIterations, centralizedCost, holds
%   mode, updateType, w1, w2, params (optional struct for gains/exponents)
%
% Notes:
%   - In 'primal' mode, hats are updated externally each ADMM iteration and
%     enter ONLY through weighted residuals in the consistency penalty terms.
%   - In 'dual' mode, hat variables and w1/w2 are unused.
% =========================================================================

clear global

% ---------------- Defaults ----------------
if ~endsWith(partitionedDataFile, '.mat')
    partitionedDataFile = [partitionedDataFile, '.mat'];
end
if nargin < 5 || isempty(centralizedCost), centralizedCost = 0; end
if nargin < 6 || isempty(holds), holds = 1e8; end
if nargin < 7 || isempty(mode), mode = 'none'; end
if nargin < 8 || isempty(updateType), updateType = ''; end
if nargin < 9 || isempty(w1), w1 = 1; end
if nargin < 10 || isempty(w2), w2 = 0; end
if nargin < 11 || isempty(params), params = struct(); end

mode = lower(string(mode));
updateType = lower(string(updateType));solverTimeTotal   = 0;

% --- Default primal gains/exponents (can override via params.*) ---
if ~isfield(params,'alpha1'), params.alpha1 = -1; end           % asy
if ~isfield(params,'alpha2'), params.alpha2 = -0.1; end           % finite/fixed
if ~isfield(params,'alpha3'), params.alpha3 = -0.05; end           % fixed
if ~isfield(params,'mu_p'),   params.mu_p   = 0.5; end            % 0<mu<1
if ~isfield(params,'eta_p'),  params.eta_p  = 2.1; end            % eta>1

% --- Default dual gains/exponents (can override via params.*) ---
if ~isfield(params,'gamma1'), params.gamma1 = 1.0; end           % finite/fixed
if ~isfield(params,'gamma2'), params.gamma2 = 1.0; end           % fixed
if ~isfield(params,'mu_d'),   params.mu_d   = 0.1; end            % 0<mu<1
if ~isfield(params,'eta_d'),  params.eta_d  = -0.1; end           % -1<eta<0
if ~isfield(params,'eps_d'),  params.eps_d  = 1e-12; end           % denom eps

% global switches for callbacks (acopf.m)
global admm_mode admm_updateType admm_w1 admm_w2 admm_params
admm_mode = mode;
admm_updateType = updateType;
admm_w1 = w1;
admm_w2 = w2;
admm_params = params;

% ---------------- Load partitioned data ----------------
global loadedData regnum
loadedData = load(partitionedDataFile);
numRegions = loadedData.num_regions;

% --- ADMM globals ---
global rho_theta rho_V PuScale landa lambdav
global V_scale theta_scale

% ---------------- Augment each region with interregional tie-lines ----------------
for regionIdx = 1:numRegions
    regionID   = ['R' num2str(regionIdx)];
    regionName = ['mpc_region', regionID];

    tielinesName = ['interregional_tielines', regionID];
    baseBranch   = loadedData.(regionName);

    dataMatrix  = table2cell(loadedData.interregional_tielines_total);
    row_indices = strcmp(dataMatrix(:, 1), regionName) | strcmp(dataMatrix(:, 2), regionName);
    filteredData = dataMatrix(row_indices, :);

    nbus = size(loadedData.(regionName).bus,1);
    checkmattrix = zeros(size(filteredData,1),1);

    for i = 1:size(filteredData,1)
        region1_str = filteredData{i, 1};
        region2_str = filteredData{i, 2};

        region1_num = sscanf(region1_str, 'mpc_regionR%d');
        region2_num = sscanf(region2_str, 'mpc_regionR%d');

        if region1_num == regionIdx
            checkmattrix(i)=1;
        else
            checkmattrix(i)=-1;
        end
    end

    output_rows = cell2mat(filteredData(:, 3:end));
    for i = 1:size(filteredData,1)
        if checkmattrix(i)==1
            output_rows(i,2)=nbus+i;
        else
            output_rows(i,1)=nbus+i;
        end
    end

    for i = 1:size(output_rows, 1)
        filteredData{i, 3} = output_rows(i, :);
    end
    loadedData.(tielinesName) = filteredData;

    % Construct struct dynamically (mpcRk)
    mpc_struct = struct;
    mpc_struct.version = '2';
    if isfield(baseBranch, 'baseMVA'), mpc_struct.baseMVA = baseBranch.baseMVA; else, mpc_struct.baseMVA = 100; end
    if isfield(baseBranch, 'bus'),     mpc_struct.bus = baseBranch.bus; end
    if isfield(baseBranch, 'gen'),     mpc_struct.gen = baseBranch.gen; end
    if isfield(baseBranch, 'branch'),  mpc_struct.branch = baseBranch.branch; end
    if isfield(baseBranch, 'gencost'), mpc_struct.gencost = baseBranch.gencost; end

    structName = ['mpc', regionID];
    assignin('base', structName, mpc_struct);

    % Compute interregional Y pieces (G,B) for tie-lines
    Gname      = ['G', regionID];
    Bname      = ['B', regionID];
    checkname  = ['check', regionID];

    LL = output_rows;
    nbranches = size(LL, 1);
    G = zeros(nbranches, 4);
    B = zeros(nbranches, 4);

    for i = 1:nbranches
        r     = LL(i, 3);
        x     = LL(i, 4);
        bc    = LL(i, 5);
        ratio = LL(i, 9);
        if ratio == 0, ratio = 1; end
        angle     = LL(i, 10);
        angle_rad = pi * angle / 180;

        invratio2 = 1 / ratio^2;
        multtf    = 1 / (ratio * exp(1j * angle_rad));
        multft    = 1 / (ratio * exp(-1j * angle_rad));
        z         = r + 1j * x;
        y         = 1 / z;

        Yff = (y + bc / 2 * 1j) * invratio2;
        Yft = -y * multft;
        Ytf = -y * multtf;
        Ytt = y + bc / 2 * 1j;

        G(i, 1) = real(Yff);  B(i, 1) = imag(Yff);
        G(i, 2) = real(Yft);  B(i, 2) = imag(Yft);
        G(i, 3) = real(Ytf);  B(i, 3) = imag(Ytf);
        G(i, 4) = real(Ytt);  B(i, 4) = imag(Ytt);
    end

    check = double(LL(:, 1) > LL(:, 2));
    loadedData.(Gname)     = G;
    loadedData.(Bname)     = B;
    loadedData.(checkname) = check;
end

% ---------------- Initialize region variables ----------------
totalInterfaces = 0;
for regionIdx = 1:numRegions
    regionID   = ['R' num2str(regionIdx)];
    regionName = ['mpc_region', regionID];
    regionData = loadedData.(regionName);
    tielinesName = ['interregional_tielines', regionID];
    tielinesData = cell2table(loadedData.(tielinesName));
    tieLines     = tielinesData.Var3;

    nb  = size(regionData.bus, 1);
    nit = size(tielinesData, 1);

    phase_angles = zeros(2 * nit, 1);
    voltage_mag  = zeros(2 * nit, 1);

    % NEW: hats (same layout as phase_angles/voltage_mag)
    hat_phase_angles = ones(2 * nit, 1);
    hat_voltage_mag  = ones(2 * nit, 1);

    Yon = [];
    con = [];

    regionBusIDs = cell(size(tieLines, 1), 1);
    for j = 1:size(tieLines, 1)
        if tieLines(j, 1) <= nb
            regionBusIDs{j} = tieLines(j, 1);
        elseif tieLines(j, 2) <= nb
            regionBusIDs{j} = tieLines(j, 2);
        end
    end
    regionBusIDs = regionBusIDs(~cellfun('isempty', regionBusIDs));

    loadedData.(['nb', regionID])             = nb;
    loadedData.(['nit', regionID])            = nit;
    loadedData.(['con', regionID])            = con;
    loadedData.(['Yon', regionID])            = Yon;
    loadedData.(['phase_angles_', regionID])  = phase_angles;
    loadedData.(['voltage_mag_', regionID])   = voltage_mag;

    % NEW hats storage
    loadedData.(['hat_phase_angles_', regionID]) = hat_phase_angles;
    loadedData.(['hat_voltage_mag_',  regionID]) = hat_voltage_mag;

    loadedData.(['region_', regionID])        = regionBusIDs;

    totalInterfaces = totalInterfaces + nit;
end

nit = totalInterfaces / 2;
loadedData.nit = nit;

% ---------------- ADMM parameters & scaling ----------------
baseMV  = loadedData.mpc_regionR1.baseMVA;
mpc0    = loadcase(loadedData.filename);
nb_full = size(mpc0.bus, 1); %#ok<NASGU>

rho_theta = rho;
rho_V     = rho;

x     = mpc0.branch(:, 4);
rateA = mpc0.branch(:, 6);
Vmin  = min(mpc0.bus(:,13));
Vmax  = max(mpc0.bus(:,12));

V_scale = (Vmax - Vmin)/2;
theta_max_approx = abs(x) .* (rateA / baseMV);
valid = theta_max_approx > 0 & ~isnan(theta_max_approx) & ~isinf(theta_max_approx);
if any(valid)
    theta_scale = median(theta_max_approx(valid));
else
    theta_scale = 0.1;
end
theta_scale = 1;V_scale = 1;

if strcmpi(mode,'none')
    PuScale = 1.0;

elseif strcmpi(mode,'primal')
    if strcmpi(updateType,'asy')
        PuScale = 1.5;
    elseif strcmpi(updateType,'finite')
        PuScale = 1.9;
    elseif strcmpi(updateType,'fixed')
        PuScale = 3.5;
    else
        error('Unknown updateType for primal: %s', updateType);
    end

elseif strcmpi(mode,'dual')
    if strcmpi(updateType,'finite')
        PuScale = 5;
    elseif strcmpi(updateType,'fixed')
        PuScale = 4.5;
    else
        error('Unknown updateType for dual: %s', updateType);
    end

else
    error('Unknown mode: %s', mode);
end

landa      = zeros(2 * nit, 1);
lambdav    = zeros(2 * nit, 1);

errorValues   = zeros(2 * nit, 1);
errorValuesv  = zeros(2 * nit, 1);
errorLog      = zeros(maxIterations, 1);
optimalityGap = zeros(maxIterations, 1);

r_theta_vec     = zeros(2 * nit, 1);
r_theta_vec_old = zeros(2 * nit, 1);
r_V_vec         = zeros(2 * nit, 1);
r_V_vec_old     = zeros(2 * nit, 1);

rho_theta_log = zeros(maxIterations,1);
rho_V_log     = zeros(maxIterations,1);

% Adaptive rho params (kept as-is for 'none'/'primal'; disabled in 'dual')
mu_theta  = max(1,15/((numRegions^2)*0.8));
mu_V      = max(1,10/((numRegions^2)*0.8));
tau_incr  = 1.5;
tau_decr  = 1/1.5;
eta_dualstep = 1; % classical dual step multiplier (your eta)
Flag_rho_theta=0; Flag_rho_v=0;
iter_count_theta=0; iter_count_v=0;
max_count = 1;

% ---------------- Prepare interregional interface mapping ----------------
tieLineTable = loadedData.interregional_tielines_total;
tieLines     = table2cell(tieLineTable);

for i = 1:nit
    r1 = tieLines{i, 1};
    r2 = tieLines{i, 2};
    r1_num = sscanf(r1, 'mpc_regionR%d');
    r2_num = sscanf(r2, 'mpc_regionR%d');
    if r1_num > r2_num
        [tieLines{i, 1}, tieLines{i, 2}] = deal(r2, r1);
    end
end
tieLineTable = cell2table(tieLines);

interfaceCounter = zeros(numRegions, 1);
for regionIdx = 1:numRegions
    regionID = ['R', num2str(regionIdx)];
    loadedData.(['con', regionID]) = [];
    loadedData.(['Yon', regionID]) = [];
end

for i = 1:nit
    region1_str = tieLines{i, 1};
    region2_str = tieLines{i, 2};
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

    reg1_id = ['R', num2str(region1_num)];
    reg2_id = ['R', num2str(region2_num)];

    loadedData.(['con', reg1_id]) = [loadedData.(['con', reg1_id]); b];
    loadedData.(['con', reg2_id]) = [loadedData.(['con', reg2_id]); a];
    loadedData.(['Yon', reg1_id]) = [loadedData.(['Yon', reg1_id]); i];
    loadedData.(['Yon', reg2_id]) = [loadedData.(['Yon', reg2_id]); i];
end

% ---------------- ADMM loop ----------------
for iter = 1:maxIterations
    currentObjective = 0;
    Pg_total = 0;

    % ===== PRIMAL MODE: update hats BEFORE local OPF solves =====
    if mode == "primal"
        for r = 1:numRegions
            regnum = r;
            regTag = ['R' num2str(r)];
            nitr   = loadedData.(['nit' regTag]);
            con    = loadedData.(['con' regTag]);
            Yon    = loadedData.(['Yon' regTag]);

            hat_phase = loadedData.(['hat_phase_angles_' regTag]);
            hat_mag   = loadedData.(['hat_voltage_mag_'  regTag]);

            % We update ONLY the COPY positions (second half): (i+nitr)
            for i = 1:nitr
                % Determine remote region for this interface i
                % Use the same table used in callbacks for robust mapping
                tieKey = ['interregional_tielines' regTag];
                dataMatrix = loadedData.(tieKey);

                % Sort numeric region order (same as callbacks)
                for q = 1:size(dataMatrix,1)
                    rr1 = dataMatrix{q,1}; rr2 = dataMatrix{q,2};
                    rr1n = sscanf(rr1,'mpc_regionR%d'); rr2n = sscanf(rr2,'mpc_regionR%d');
                    if rr1n > rr2n
                        dataMatrix{q,1}=rr2; dataMatrix{q,2}=rr1;
                    end
                end
                Tablek = cell2table(dataMatrix);

                str1 = char(Tablek{i,1});
                str2 = char(Tablek{i,2});
                region1_num = sscanf(str1,'mpc_regionR%d');
                region2_num = sscanf(str2,'mpc_regionR%d');
                localRegion = r;
                if localRegion == region1_num
                    remoteRegion = region2_num;
                else
                    remoteRegion = region1_num;
                end
                remoteTag = ['R' num2str(remoteRegion)];

                remotePhase = loadedData.(['phase_angles_' remoteTag]);
                remoteMag   = loadedData.(['voltage_mag_'  remoteTag]);

                % remote owned voltage for this interface is at index con(i)
                ref_theta = remotePhase(con(i));
                ref_v     = remoteMag(con(i));

                % current hat (copy slot)
                idxHat = i + nitr;
                e_theta = hat_phase(idxHat) - ref_theta;
                e_v     = hat_mag(idxHat)   - ref_v;

                switch updateType
                    case "asy"
                        hat_phase(idxHat) = hat_phase(idxHat) + params.alpha1 * e_theta;
                        hat_mag(idxHat)   = hat_mag(idxHat)   + params.alpha1 * e_v;

                    case "finite"
                        mu_p = params.mu_p;
                        hat_phase(idxHat) = hat_phase(idxHat) + params.alpha2 * sign(e_theta) * (abs(e_theta)^mu_p);
                        hat_mag(idxHat)   = hat_mag(idxHat)   + params.alpha2 * sign(e_v)     * (abs(e_v)^mu_p);

                    case "fixed"
                        mu_p  = params.mu_p;
                        eta_p = params.eta_p;
                        hat_phase(idxHat) = hat_phase(idxHat) ...
                            + params.alpha3 * sign(e_theta) * (abs(e_theta)^mu_p) ...
                            + params.alpha3 * sign(e_theta) * (abs(e_theta)^eta_p);
                        hat_mag(idxHat)   = hat_mag(idxHat) ...
                            + params.alpha3 * sign(e_v)     * (abs(e_v)^mu_p) ...
                            + params.alpha3 * sign(e_v)     * (abs(e_v)^eta_p);

                    otherwise
                        error("In primal mode, updateType must be 'asy', 'finite', or 'fixed'.");
                end
            end

            loadedData.(['hat_phase_angles_' regTag]) = hat_phase;
            loadedData.(['hat_voltage_mag_'  regTag]) = hat_mag;
        end
    end

    % ===== Local ACOPF solves =====
    for r = 1:numRegions
        reg    = ['R' num2str(r)];
        regnum = r;

        nb      = loadedData.(['nb' reg]);
        nitall  = loadedData.(['nit' reg]);
        regionName = loadedData.(['region_' reg]);

        [Optimized_Solution, gen_coeff, flag, obj, stat,ytyt] = acopf(reg); %#ok<ASGLU>
        solverTimeTotal = solverTimeTotal + ytyt;
        nvars = numel(Optimized_Solution) - 2 * nitall;

        phase_angles = zeros(2 * nitall, 1);
        voltage_mag  = zeros(2 * nitall, 1);

        for i = 1:nitall
            phase_angles(i, 1)        = Optimized_Solution(regionName{i} + nb);
            phase_angles(i+nitall, 1) = Optimized_Solution((2 * i) + nvars);

            voltage_mag(i, 1)         = Optimized_Solution(regionName{i});
            voltage_mag(i+nitall, 1)  = Optimized_Solution(((2 * i) - 1) + nvars);
        end

        loadedData.(['phase_angles_', reg]) = phase_angles;
        loadedData.(['voltage_mag_',  reg]) = voltage_mag;

        if flag ~= 0
            Pg = zeros(size(gen_coeff, 1), 1);
            for i = 1:numel(Pg)
                Pg(i) = Optimized_Solution(2 * nb + i);
            end
        else
            Pg = zeros(size(gen_coeff,1),1);
        end

        currentObjective = currentObjective + value(obj);
        Pg_total         = Pg_total + sum(value(Pg));

        % ===== Initialize hats once (iter==1) to match copies (safe start) =====
        if (mode == "primal") && (iter == 1)
            regTag = reg;
            nitr = loadedData.(['nit' regTag]);
            hat_phase = loadedData.(['hat_phase_angles_' regTag]);
            hat_mag   = loadedData.(['hat_voltage_mag_'  regTag]);

            % copy positions are second half: i+nitr
            hat_phase((1:nitr)+nitr) = phase_angles((1:nitr)+nitr);
            hat_mag((1:nitr)+nitr)   = voltage_mag((1:nitr)+nitr);

            loadedData.(['hat_phase_angles_' regTag]) = hat_phase;
            loadedData.(['hat_voltage_mag_'  regTag]) = hat_mag;
        end
    end

    % ===== Dual update step =====
    interfaceCounter = zeros(numRegions, 1);

    for i = 1:nit
        regionName1 = char(tieLines{i, 1});
        regionName2 = char(tieLines{i, 2});

        regionIdx1 = sscanf(regionName1, 'mpc_regionR%d');
        regionIdx2 = sscanf(regionName2, 'mpc_regionR%d');

        for regionIdx = 1:numRegions
            if regionIdx == regionIdx1
                interfaceCounter(regionIdx) = interfaceCounter(regionIdx) + 1;
                interfacePos1 = interfaceCounter(regionIdx);
            elseif regionIdx == regionIdx2
                interfaceCounter(regionIdx) = interfaceCounter(regionIdx) + 1;
                interfacePos2 = interfaceCounter(regionIdx);
            end
        end

        reg1 = ['R' num2str(regionIdx1)];
        reg2 = ['R' num2str(regionIdx2)];

        anglesReg1 = loadedData.(['phase_angles_' reg1]);
        anglesReg2 = loadedData.(['phase_angles_' reg2]);
        magsReg1   = loadedData.(['voltage_mag_'  reg1]);
        magsReg2   = loadedData.(['voltage_mag_'  reg2]);
        nitReg1    = loadedData.(['nit' reg1]);
        nitReg2    = loadedData.(['nit' reg2]);

        angleGrad     = anglesReg1(interfacePos1)             - anglesReg2(nitReg2 + interfacePos2);
        angleGradDual = anglesReg1(nitReg1 + interfacePos1)   - anglesReg2(interfacePos2);

        magGrad       = magsReg1(interfacePos1)               - magsReg2(nitReg2 + interfacePos2);
        magGradDual   = magsReg1(nitReg1 + interfacePos1)     - magsReg2(interfacePos2);

        % normalized residuals (these are g = ∇_λ L in your paper)
        g_theta1 = angleGrad     / theta_scale;
        g_theta2 = angleGradDual / theta_scale;
        g_v1     = magGrad       / V_scale;
        g_v2     = magGradDual   / V_scale;

        r_theta_vec(i)      = g_theta1;
        r_theta_vec(i+nit)  = g_theta2;
        r_V_vec(i)          = g_v1;
        r_V_vec(i+nit)      = g_v2;

        if mode == "dual"
            % finite/fixed-time dual updates (paper)
            eps0 = params.eps_d;
            mu_d = params.mu_d;
            rhobartheta = rho_theta/(rho_theta^mu_d);
            rhobarv = rho_V/(rho_V^mu_d);

            switch updateType
                case "finite"
                    landa(i)       = landa(i)      + params.gamma1 * rhobartheta * (g_theta1 / (max(abs(g_theta1), eps0)^mu_d));
                    landa(i+nit)   = landa(i+nit)  + params.gamma1 * rhobartheta * (g_theta2 / (max(abs(g_theta2), eps0)^mu_d));
                    lambdav(i)     = lambdav(i)    + params.gamma1 * rhobarv * (g_v1     / (max(abs(g_v1),     eps0)^mu_d));
                    lambdav(i+nit) = lambdav(i+nit)+ params.gamma1 * rhobarv * (g_v2     / (max(abs(g_v2),     eps0)^mu_d));

                case "fixed"
                    eta_d = params.eta_d; % -1 < eta < 0
                    landa(i)       = landa(i)      + params.gamma2 * rhobartheta * (g_theta1 / (max(abs(g_theta1), eps0)^mu_d)) ...
                                               + params.gamma2 * rhobartheta * (g_theta1 / (max(abs(g_theta1), eps0)^eta_d));
                    landa(i+nit)   = landa(i+nit)  + params.gamma2 * rhobartheta * (g_theta2 / (max(abs(g_theta2), eps0)^mu_d)) ...
                                               + params.gamma2 * rhobartheta * (g_theta2 / (max(abs(g_theta2), eps0)^eta_d));
                    lambdav(i)     = lambdav(i)    + params.gamma2 * rhobarv * (g_v1     / (max(abs(g_v1),     eps0)^mu_d)) ...
                                               + params.gamma2 * rhobarv * (g_v1     / (max(abs(g_v1),     eps0)^eta_d));
                    lambdav(i+nit) = lambdav(i+nit)+ params.gamma2 * rhobarv * (g_v2     / (max(abs(g_v2),     eps0)^mu_d)) ...
                                               + params.gamma2 * rhobarv * (g_v2     / (max(abs(g_v2),     eps0)^eta_d));
                otherwise
                    error("In dual mode, updateType must be 'finite' or 'fixed'.");
            end

        else
            % classical dual update (your current)
            landa(i)        = landa(i)        + eta_dualstep * rho_theta * g_theta1;
            landa(i + nit)  = landa(i + nit)  + eta_dualstep * rho_theta * g_theta2;
            lambdav(i)      = lambdav(i)      + eta_dualstep * rho_V     * g_v1;
            lambdav(i + nit)= lambdav(i + nit)+ eta_dualstep * rho_V     * g_v2;
        end

        errorValues(i)        = abs(g_theta1);
        errorValues(i + nit)  = abs(g_theta2);
        errorValuesv(i)       = abs(g_v1);
        errorValuesv(i + nit) = abs(g_v2);
    end

    % ===== Logging & convergence =====
    maxError = max([max(errorValues), max(errorValuesv)]);
    errorLog(iter) = maxError;

    if centralizedCost ~= 0
        optimalityGap(iter) = abs((currentObjective - centralizedCost) / centralizedCost) * 100;
    else
        optimalityGap(iter) = NaN;
    end

    primal_theta = max(abs(r_theta_vec));
    primal_V     = max(abs(r_V_vec));
    dual_theta   = norm(r_theta_vec - r_theta_vec_old, inf);
    dual_V       = norm(r_V_vec     - r_V_vec_old,     inf);

    % Adaptive rho only for 'none'/'primal' (disable in 'dual')
    
        if (Flag_rho_theta==0 && rho_theta<holds)
            if primal_theta > mu_theta * dual_theta
                rho_theta = rho_theta * tau_incr; Flag_rho_theta=1;
            elseif dual_theta > mu_theta * primal_theta
                rho_theta = rho_theta / tau_decr; Flag_rho_theta=1;
            end
        else
            iter_count_theta=iter_count_theta+1;
        end
        if iter_count_theta==max_count, iter_count_theta=1; Flag_rho_theta=0; end

        if (Flag_rho_v==0 && rho_V<holds)
            if primal_V > mu_V * dual_V
                rho_V = rho_V * tau_incr; Flag_rho_v=1;
            elseif dual_V > mu_V * primal_V
                rho_V = rho_V / tau_decr; Flag_rho_v=1;
            end
        else
            iter_count_v=iter_count_v+1;
        end
        if iter_count_v==max_count, iter_count_v=1; Flag_rho_v=0; end
    

    r_theta_vec_old = r_theta_vec;
    r_V_vec_old     = r_V_vec;
    rho_theta_log(iter) = rho_theta;
    rho_V_log(iter)     = rho_V;

    if maxError < residualThreshold
        fprintf('Convergence achieved at iteration %d\n', iter);
        if centralizedCost ~= 0
            fprintf('Iter %3d | Error: %.6e | Obj: %.4f | Gap: %.3f %% | rho_theta=%.3g | rho_V=%.3g | mode=%s | type=%s\n', ...
                iter, maxError, currentObjective, optimalityGap(iter), rho_theta, rho_V, mode, updateType);
        else
            fprintf('Iter %3d | Error: %.6e | Obj: %.4f | rho_theta=%.3g | rho_V=%.3g | mode=%s | type=%s\n', ...
                iter, maxError, currentObjective, rho_theta, rho_V, mode, updateType);
        end
        break;
    end

    if centralizedCost ~= 0
        fprintf('Iter %3d | Error: %.6e | Obj: %.4f | Gap: %.3f %% | rho_theta=%.3g | rho_V=%.3g | mode=%s | type=%s\n', ...
            iter, maxError, currentObjective, optimalityGap(iter), rho_theta, rho_V, mode, updateType);
    else
        fprintf('Iter %3d | Obj = %.6f | MaxError = %.6f | rho_theta=%.3g | rho_V=%.3g | mode=%s | type=%s\n', ...
            iter, currentObjective, maxError, rho_theta, rho_V, mode, updateType);
    end
end

if iter >= maxIterations
    disp('Maximum iterations reached without full convergence');
end
fprintf('Final Objective after %d iterations = %.6f\n', iter, currentObjective);

% ---------------- Demand/Generation summary ----------------
mpc = loadcase(loadedData.filename);
Pd  = mpc.bus(:, 3);
total_Pd = sum(Pd);
total_Pg = (mpc.baseMVA).*Pg_total;

fprintf('Total Real Power Demand (Pd): %.2f MW\n', total_Pd);
fprintf('Total Real Power Generation (Pg): %.2f MW\n', total_Pg);

% Report total wall time spent solving all local DCOPF subproblems
fprintf('Total sum of ACOPF subproblem solver time: %.4f seconds\n', ...
    solverTimeTotal);

% ---------------- Pack outputs ----------------
out = struct();
out.partitionedDataFile = partitionedDataFile;
out.rho_input           = rho;
out.rho_theta_final     = rho_theta;
out.rho_V_final         = rho_V;
out.residualThreshold   = residualThreshold;
out.maxIterations       = maxIterations;
out.iterationsRun       = iter;
out.errorLog            = errorLog(1:iter);
out.optimalityGap       = optimalityGap(1:iter);
out.rho_theta_log       = rho_theta_log(1:iter);
out.rho_V_log           = rho_V_log(1:iter);
out.landa               = landa;
out.lambdav             = lambdav;
out.total_Pd            = total_Pd;
out.total_Pg            = total_Pg;
out.currentObjective    = currentObjective;
out.loadedData          = loadedData;

% Echo mode config
out.mode        = char(mode);
out.updateType  = char(updateType);
out.w1          = w1;
out.w2          = w2;
out.params      = params;

end
