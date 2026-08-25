% -------------------------------------------------------------------------
function g = gradient(x, auxdata)
    mpc = auxdata{1};
    nbuses = size(mpc.bus, 1);
    ngens  = size(mpc.gen, 1);
    baseMVA = mpc.baseMVA;
    Pg = 2 * nbuses + (1:ngens);

    global regnum loadedData PuScale landa lambdav
    global V_scale theta_scale
    global rho_theta rho_V
    global admm_mode admm_w1 admm_w2

    regionTag = ['R' num2str(regnum)];
    nitr   = loadedData.(['nit', regionTag]);
    nitall = loadedData.('nit');

    g = zeros(2 * nbuses + 2 * ngens + 2 * nitr, 1);

    % generation cost gradient (unchanged)
    if ngens~=0
        actgen = mpc.gen(:, 8);
        if mpc.gencost(1, 4) == 3
            g(Pg) = actgen .* (2 * baseMVA^2 * mpc.gencost(:, 5) .* x(Pg) + baseMVA * mpc.gencost(:, 6));
        elseif mpc.gencost(1, 4) == 2
            g(Pg) = actgen .* (baseMVA * mpc.gencost(:, 5));
        end
    end

    tieTableKey = ['interregional_tielines', regionTag];
    dataMatrix = loadedData.(tieTableKey);

    for i = 1:size(dataMatrix, 1)
        r1 = dataMatrix{i, 1};
        r2 = dataMatrix{i, 2};
        r1_num = sscanf(r1, 'mpc_regionR%d');
        r2_num = sscanf(r2, 'mpc_regionR%d');
        if r1_num > r2_num
            dataMatrix{i, 1} = r2;
            dataMatrix{i, 2} = r1;
        end
    end
    Tablek = cell2table(dataMatrix);

    con   = loadedData.(['con', regionTag]);
    Yon   = loadedData.(['Yon', regionTag]);
    region_buses = loadedData.(['region_', regionTag]);

    hat_phase = loadedData.(['hat_phase_angles_' regionTag]);
    hat_mag   = loadedData.(['hat_voltage_mag_'  regionTag]);

    if admm_mode == "primal"
        denom = (admm_w1 + admm_w2); if denom <= 0, denom = 1; end
        cw = admm_w1 / denom;   % IMPORTANT: derivative scaling for copy vars
    else
        cw = 1;
    end

    for i = 1:nitr
        str1 = char(Tablek{i, 1});
        str2 = char(Tablek{i, 2});
        region1_num = sscanf(str1, 'mpc_regionR%d');
        region2_num = sscanf(str2, 'mpc_regionR%d');
        localRegion_num = regnum;

        if localRegion_num == region1_num
            remoteRegion_num = region2_num;
        else
            remoteRegion_num = region1_num;
        end
        remoteTag = ['R' num2str(remoteRegion_num)];

        remotePhase = loadedData.(['phase_angles_', remoteTag]);
        remoteMag   = loadedData.(['voltage_mag_',  remoteTag]);
        nit_remote  = loadedData.(['nit', remoteTag]);

        % local vars
        owned_theta = x(region_buses{i} + nbuses);
        owned_v     = x(region_buses{i});
        copy_theta  = x(2*nbuses + 2*ngens + 2*i);
        copy_v      = x(2*nbuses + 2*ngens + 2*i - 1);

        % remote constants
        remote_owned_theta = remotePhase(con(i));
        remote_owned_v     = remoteMag(con(i));
        remote_copy_theta  = remotePhase(nit_remote + con(i));
        remote_copy_v      = remoteMag(nit_remote + con(i));

        if admm_mode == "primal"
            denom = (admm_w1 + admm_w2); if denom <= 0, denom = 1; end
            copy_theta_eff = (admm_w1 * copy_theta + admm_w2 * hat_phase(i+nitr)) / denom;
            copy_v_eff     = (admm_w1 * copy_v     + admm_w2 * hat_mag(i+nitr))   / denom;
        else
            copy_theta_eff = copy_theta;
            copy_v_eff     = copy_v;
        end

        y = Yon(i);

        % residuals (same as objective)
        if localRegion_num == region1_num
            r_theta      = (owned_theta - remote_copy_theta) / theta_scale;          % owned var
            r_theta_dual = (copy_theta_eff - remote_owned_theta) / theta_scale;      % copy var (scaled by cw)
            r_v          = (owned_v     - remote_copy_v) / V_scale;                  % owned var
            r_v_dual     = (copy_v_eff  - remote_owned_v) / V_scale;                 % copy var (scaled by cw)
        else
            r_theta      = (remote_owned_theta - copy_theta_eff) / theta_scale;      % copy var (scaled by cw)
            r_theta_dual = (remote_copy_theta  - owned_theta) / theta_scale;         % owned var
            r_v          = (remote_owned_v     - copy_v_eff) / V_scale;              % copy var (scaled by cw)
            r_v_dual     = (remote_copy_v      - owned_v) / V_scale;                 % owned var
        end

        idxOwnedV     = region_buses{i};
        idxOwnedTheta = region_buses{i} + nbuses;
        idxCopyV      = 2*nbuses + 2*ngens + 2*i - 1;
        idxCopyTheta  = 2*nbuses + 2*ngens + 2*i;

        if ngens~=0
            s = PuScale;
        else
            s = 1;
        end

        % contributions:
        % owned variables: derivative coefficient = 1
        % copy variables (if primal): derivative coefficient = cw
        if localRegion_num == region1_num
            % r_theta depends on owned_theta
            g(idxOwnedTheta) = g(idxOwnedTheta) + s * (landa(abs(y)) + rho_theta * r_theta) / theta_scale;

            % r_theta_dual depends on copy_theta_eff => derivative wrt copy_theta is cw
            g(idxCopyTheta)  = g(idxCopyTheta)  + s * cw * (landa(abs(y)+nitall) + rho_theta * r_theta_dual) / theta_scale;

            % r_v depends on owned_v
            g(idxOwnedV)     = g(idxOwnedV)     + s * (lambdav(abs(y)) + rho_V * r_v) / V_scale;

            % r_v_dual depends on copy_v_eff => derivative wrt copy_v is cw
            g(idxCopyV)      = g(idxCopyV)      + s * cw * (lambdav(abs(y)+nitall) + rho_V * r_v_dual) / V_scale;
        else
            % r_theta depends on copy_theta_eff with NEGATIVE sign; derivative wrt copy is -cw
            g(idxCopyTheta)  = g(idxCopyTheta)  + s * (-cw) * (landa(abs(y)) + rho_theta * r_theta) / theta_scale;

            % r_theta_dual depends on owned_theta with NEGATIVE sign
            g(idxOwnedTheta) = g(idxOwnedTheta) + s * (-1)  * (landa(abs(y)+nitall) + rho_theta * r_theta_dual) / theta_scale;

            % r_v depends on copy_v_eff with NEGATIVE sign
            g(idxCopyV)      = g(idxCopyV)      + s * (-cw) * (lambdav(abs(y)) + rho_V * r_v) / V_scale;

            % r_v_dual depends on owned_v with NEGATIVE sign
            g(idxOwnedV)     = g(idxOwnedV)     + s * (-1)  * (lambdav(abs(y)+nitall) + rho_V * r_v_dual) / V_scale;
        end
    end