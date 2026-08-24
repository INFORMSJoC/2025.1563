% -------------------------------------------------------------------------
function f = objective(x, auxdata)
    mpc = auxdata{1};
    nbuses = size(mpc.bus, 1);
    ngens  = size(mpc.gen, 1);
    baseMVA = mpc.baseMVA;
    Pg = 2*nbuses + (1:ngens);

    % Base generation cost (unchanged)
    if ngens~=0
        actgen = mpc.gen(:, 8);
        if mpc.gencost(1, 4) == 3
            f = sum(actgen .* (baseMVA^2 * mpc.gencost(:, 5) .* x(Pg).^2 + baseMVA * mpc.gencost(:, 6) .* x(Pg) + mpc.gencost(:, 7)));
        elseif mpc.gencost(1, 4) == 2
            f = sum(actgen .* (baseMVA * mpc.gencost(:, 5) .* x(Pg) + mpc.gencost(:, 6)));
        end
    else
        f = 0;
    end

    % Consistency penalties (modified only in primal mode for copy-residuals)
    global PuScale landa lambdav regnum loadedData
    global V_scale theta_scale
    global rho_theta rho_V
    global admm_mode admm_w1 admm_w2

    regionTag = ['R' num2str(regnum)];
    tieTableKey = ['interregional_tielines' regionTag];
    dataMatrix = loadedData.(tieTableKey);

    % numeric sort (unchanged)
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

    con   = loadedData.(['con' regionTag]);
    Yon   = loadedData.(['Yon' regionTag]);
    region_buses = loadedData.(['region_' regionTag]);
    nitr  = loadedData.(['nit' regionTag]);
    nitall= loadedData.('nit');

    % hats (only used in primal mode)
    hat_phase = loadedData.(['hat_phase_angles_' regionTag]);
    hat_mag   = loadedData.(['hat_voltage_mag_'  regionTag]);


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
        remoteRegionStr = ['R' num2str(remoteRegion_num)];

        remotePhase = loadedData.(['phase_angles_' remoteRegionStr]);
        remoteMag   = loadedData.(['voltage_mag_'  remoteRegionStr]);
        remote_nit  = loadedData.(['nit' remoteRegionStr]);

        y = Yon(i);

        % Local variables
        owned_theta = x(region_buses{i} + nbuses);
        owned_v     = x(region_buses{i});
        copy_theta  = x(2*nbuses + 2*ngens + 2*i);
        copy_v      = x(2*nbuses + 2*ngens + 2*i - 1);

        % Remote constants
        remote_owned_theta = remotePhase(con(i));
        remote_owned_v     = remoteMag(con(i));
        remote_copy_theta  = remotePhase(remote_nit + con(i));
        remote_copy_v      = remoteMag(remote_nit + con(i));

        % Weighted copy value (only in primal mode)
        if admm_mode == "primal"
            denom = (admm_w1 + admm_w2); if denom <= 0, denom = 1; end
            copy_theta_eff = (admm_w1 * copy_theta + admm_w2 * hat_phase(i+nitr)) / denom;
            copy_v_eff     = (admm_w1 * copy_v     + admm_w2 * hat_mag(i+nitr))   / denom;
        else
            copy_theta_eff = copy_theta;
            copy_v_eff     = copy_v;
        end

        % Residuals (same structure as your original, but with copy_eff when needed)
        if localRegion_num == region1_num
            r_theta      = (owned_theta - remote_copy_theta) / theta_scale;   % untouched
            r_theta_dual = (copy_theta_eff - remote_owned_theta) / theta_scale; % modified in primal
            r_v          = (owned_v     - remote_copy_v) / V_scale;           % untouched
            r_v_dual     = (copy_v_eff  - remote_owned_v) / V_scale;          % modified in primal
        else
            r_theta      = (remote_owned_theta - copy_theta_eff) / theta_scale; % modified in primal
            r_theta_dual = (remote_copy_theta  - owned_theta) / theta_scale;   % untouched
            r_v          = (remote_owned_v     - copy_v_eff) / V_scale;        % modified in primal
            r_v_dual     = (remote_copy_v      - owned_v) / V_scale;           % untouched
        end

        % Augmented Lagrangian accumulation (unchanged structure)
        if ngens~=0
            f = f + PuScale * ( ...
                  landa(abs(y))            * r_theta      + (rho_theta / 2) * r_theta^2 + ...
                  landa(abs(y) + nitall)   * r_theta_dual + (rho_theta / 2) * r_theta_dual^2 + ...
                  lambdav(abs(y))          * r_v          + (rho_V / 2)     * r_v^2 + ...
                  lambdav(abs(y) + nitall) * r_v_dual     + (rho_V / 2)     * r_v_dual^2 );
        else
            f = f + ( ...
                  landa(abs(y))            * r_theta      + (rho_theta / 2) * r_theta^2 + ...
                  landa(abs(y) + nitall)   * r_theta_dual + (rho_theta / 2) * r_theta_dual^2 + ...
                  lambdav(abs(y))          * r_v          + (rho_V / 2)     * r_v^2 + ...
                  lambdav(abs(y) + nitall) * r_v_dual     + (rho_V / 2)     * r_v_dual^2 );
        end
    end

