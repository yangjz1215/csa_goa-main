% calcFitness.m - 适应度计算函数
% 参考 Code_main0\Objective.m
% 适应度 = 满足实时约束的任务优先级加权和
% UAV_type: 0=普通, 1=增强型(eUAV)
% RRH_type: 0=普通, 1=增强型(eRRH)
function [fit, weighted_success] = calcFitness(uav_pos, User, priorities, E_remaining, E_max, k_move, g, subpop_params, N_UAV, cover_radius, RRH, capturability, N_RRH, RRH_type, UAV_type, params)
    uav_pos = reshape(uav_pos, N_UAV, 2);
    N_User = size(User, 1);
    
    H = 50;
    omega0 = 1e6;
    alpha0 = 1.42e-4;
    sigma2 = 3.98e-12;
    f_eRRH = 2e9;
    f_eUAV = 2e9;
    f_BBU = 4e9;
    Ptx = 1;
    PtxR = 10;
    PtxU = 10;
    PtxEU = 5;
    
    if nargin >= 16 && isfield(params, 'f_eRRH')
        f_eRRH = params.f_eRRH;
        f_eUAV = params.f_eUAV;
        f_BBU = params.f_BBU;
        Ptx = params.Ptx;
        PtxR = params.PtxR;
        PtxU = params.PtxU;
        PtxEU = params.PtxEU;
        alpha0 = params.alpha0;
        sigma2 = params.sigma2;
    end
    
    % 计算当前部署方案的初始到达能耗（假设从地图中心 [500, 500] 起飞）
    center_point = repmat([500, 500], N_UAV, 1);
    fly_distances = sqrt(sum((uav_pos - center_point).^2, 2));
    fly_energy_cost = k_move * fly_distances;
    actual_E_remaining = max(0, E_max(:) - fly_energy_cost(:));
    
    if nargin >= 13 && ~isempty(N_RRH)
        Connect = UserConnect_EnergyAware(User, N_User, RRH, N_RRH, RRH_type, uav_pos, N_UAV, UAV_type, priorities, actual_E_remaining, E_max, params);
    else
        Connect = zeros(1, N_User);
        D_RRH = 150;
        D_UAV = 150;
        for i = 1:N_User
            dists_rrh = sqrt(sum((RRH - User(i,:)).^2, 2));
            dists_uav = sqrt(sum((uav_pos - User(i,:)).^2, 2) + H^2);
            if min(dists_rrh) <= D_RRH
                Connect(i) = find(dists_rrh <= D_RRH, 1);
            elseif min(dists_uav) <= D_UAV
                Connect(i) = N_RRH + find(dists_uav <= D_UAV, 1);
            end
        end
    end
    
    if nargin >= 16 && ~isempty(params)
        off_decision = UserOffloading_EnergyAware(User, N_User, N_RRH, RRH, RRH_type, uav_pos, N_UAV, UAV_type, ...
            params.D, params.C, priorities, Connect, actual_E_remaining, E_max, params);
    else
        off_decision = zeros(1, N_User);
        off_decision(Connect > 0) = Connect(Connect > 0);
    end
    
    constraint_penalty = 0;
    
    a = zeros(1, N_User);
    latency = zeros(1, N_User);
    
    for i = 1:N_User
        if Connect(i) ~= 0
            index = off_decision(i);
            if index ~= 0
                if index <= N_RRH
                    h0 = alpha0/sum((RRH(index,:)-User(i,:)).^2);
                    r = omega0*log2(1+(Ptx*h0)/sigma2);
                    
                    if isvector(params.D) && length(params.D) > 1
                        D_user = params.D(i);
                    else
                        D_user = params.D;
                    end
                    t_tx = D_user/r;
                    
                    if isvector(params.C) && length(params.C) > 1
                        C_user = params.C(i);
                    else
                        C_user = params.C;
                    end
                    t_ex = C_user/f_eRRH;
                    latency(i) = t_ex + t_tx;
                    
                    if isvector(params.DT) && length(params.DT) > 1
                        DT_user = params.DT(i);
                    else
                        DT_user = params.DT;
                    end
                    
                    if latency(i) <= DT_user
                        a(i) = 1;
                    end
                    
                elseif index > N_RRH && index <= N_RRH + N_UAV
                    uav_idx = index - N_RRH;
                    h0 = alpha0/sum((uav_pos(uav_idx,:)-User(i,:)).^2 + H^2);
                    r = omega0*log2(1+(Ptx*h0)/sigma2);
                    
                    if isvector(params.D) && length(params.D) > 1
                        D_user = params.D(i);
                    else
                        D_user = params.D;
                    end
                    t_tx = D_user/r;
                    
                    if isvector(params.C) && length(params.C) > 1
                        C_user = params.C(i);
                    else
                        C_user = params.C;
                    end
                    t_ex = C_user/f_eUAV;
                    latency(i) = t_ex + t_tx;
                    
                    if isvector(params.DT) && length(params.DT) > 1
                        DT_user = params.DT(i);
                    else
                        DT_user = params.DT;
                    end
                    
                    if latency(i) <= DT_user
                        a(i) = 1;
                    end
                    
                else
                    cindex = Connect(i);
                    if cindex <= N_RRH
                        h0 = alpha0/sum((RRH(cindex,:)-User(i,:)).^2);
                        r1 = omega0*log2(1+(Ptx*h0)/sigma2);
                        
                        if isvector(params.D) && length(params.D) > 1
                            D_user = params.D(i);
                        else
                            D_user = params.D;
                        end
                        t_tx1 = D_user/r1;
                        
                        h01 = alpha0/sum((RRH(cindex,:)).^2);
                        r2 = omega0*log2(1+(PtxR*h01)/sigma2);
                        t_tx2 = D_user/r2;
                        
                        if isvector(params.C) && length(params.C) > 1
                            C_user = params.C(i);
                        else
                            C_user = params.C;
                        end
                        t_ex = C_user/f_BBU;
                        latency(i) = t_ex + t_tx1 + t_tx2;
                        
                        if isvector(params.DT) && length(params.DT) > 1
                            DT_user = params.DT(i);
                        else
                            DT_user = params.DT;
                        end
                        
                        if latency(i) <= DT_user
                            a(i) = 1;
                        end
                    else
                        uav_idx = cindex - N_RRH;
                        if UAV_type(uav_idx) == 1
                            PtxUU = PtxEU;
                        else
                            PtxUU = PtxU;
                        end
                        
                        h0 = alpha0/sum((uav_pos(uav_idx,:)-User(i,:)).^2 + H^2);
                        r1 = omega0*log2(1+(Ptx*h0)/sigma2);
                        
                        if isvector(params.D) && length(params.D) > 1
                            D_user = params.D(i);
                        else
                            D_user = params.D;
                        end
                        t_tx1 = D_user/r1;
                        
                        h01 = alpha0/sum((uav_pos(uav_idx,:)).^2 + H^2);
                        r2 = omega0*log2(1+(PtxUU*h01)/sigma2);
                        t_tx2 = D_user/r2;
                        
                        if isvector(params.C) && length(params.C) > 1
                            C_user = params.C(i);
                        else
                            C_user = params.C;
                        end
                        t_ex = C_user/f_BBU;
                        latency(i) = t_ex + t_tx1 + t_tx2;
                        
                        if isvector(params.DT) && length(params.DT) > 1
                            DT_user = params.DT(i);
                        else
                            DT_user = params.DT;
                        end
                        
                        if latency(i) <= DT_user
                            a(i) = 1;
                        end
                    end
                end
            end
        end
    end
    
    weighted_success = sum(priorities(a == 1));
    base_fitness = weighted_success;
    
    priorities_sum = sum(priorities);
    if priorities_sum == 0
        priorities_sum = 1;
    end
    
    high_mask = priorities >= 3;
    low_mask = priorities < 3;
    if any(high_mask)
        Succ_high = mean(a(high_mask) == 1);
    else
        Succ_high = mean(a == 1);
    end
    if any(low_mask)
        Succ_low = mean(a(low_mask) == 1);
    else
        Succ_low = mean(a == 1);
    end
    Succ_total = mean(a == 1);
    
    if g == 3
        % 使用 actual_E_remaining（已考虑飞行能耗），计算归一化能耗
        total_energy_consumed = sum(E_max - actual_E_remaining);
        normalized_energy = total_energy_consumed / (N_UAV * E_max);
    else
        normalized_energy = 0;
    end
    
    uav_dispersion = 0;
    if N_UAV > 1
        min_distances = zeros(N_UAV, 1);
        for i = 1:N_UAV
            dists = sqrt(sum((uav_pos - repmat(uav_pos(i,:), N_UAV, 1)).^2, 2));
            dists(i) = Inf;
            min_distances(i) = min(dists);
        end
        uav_dispersion = 1.5 * mean(min_distances) / (params.D_UU * 2);
        uav_dispersion = min(1, uav_dispersion);
    end

    distance_penalty = 0;
    if N_UAV > 1 && isfield(params, 'D_UU')
        for i = 1:N_UAV
            for j = i+1:N_UAV
                dist = norm(uav_pos(i,:) - uav_pos(j,:));
                if dist < params.D_UU
                    distance_penalty = distance_penalty + 30 * (params.D_UU - dist) / params.D_UU;
                end
            end
        end
    end

    base_w = 0.95;
    
    if g == 1
        fit = base_w * base_fitness ...
            + 0.03 * priorities_sum * Succ_low ...
            + 0.02 * priorities_sum * uav_dispersion ...
            - distance_penalty;
    elseif g == 2
        fit = base_w * base_fitness ...
            + 0.05 * priorities_sum * uav_dispersion ...
            - distance_penalty;
    else
        fit = base_w * base_fitness ...
            + 0.04 * priorities_sum * (1 - normalized_energy) ...
            + 0.01 * priorities_sum * uav_dispersion ...
            - distance_penalty;
    end
    fit = max(0, double(fit(1)));
end
