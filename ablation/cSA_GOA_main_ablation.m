function [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = cSA_GOA_main_ablation(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities, variant)
    if nargin < 13
        variant = 'proposed';
    end

    params.variant = variant;

    switch variant
        case 'proposed'
            params.enable_early_stop = true;
            params.enable_targeted_levy = true;
            params.enable_multi_subpop = true;
            params.enable_isolation_penalty = true;
            params.enable_smart_stop = true;
            params.enable_bayesian = true;
            params.enable_adaptive = true;

        case 'no_levy'
            params.enable_early_stop = true;
            params.enable_targeted_levy = false;
            params.enable_multi_subpop = true;
            params.enable_isolation_penalty = true;
            params.enable_smart_stop = true;
            params.enable_bayesian = true;
            params.enable_adaptive = true;

        case 'no_subpop'
            params.enable_early_stop = true;
            params.enable_targeted_levy = true;
            params.enable_multi_subpop = false;
            params.enable_isolation_penalty = true;
            params.enable_smart_stop = true;
            params.enable_bayesian = true;
            params.enable_adaptive = true;

        case 'no_penalty'
            params.enable_early_stop = true;
            params.enable_targeted_levy = true;
            params.enable_multi_subpop = true;
            params.enable_isolation_penalty = false;
            params.enable_smart_stop = true;
            params.enable_bayesian = true;
            params.enable_adaptive = true;

        case 'no_stop'
            params.enable_early_stop = false;
            params.enable_targeted_levy = true;
            params.enable_multi_subpop = true;
            params.enable_isolation_penalty = true;
            params.enable_smart_stop = false;
            params.enable_bayesian = true;
            params.enable_adaptive = true;

        case 'no_bayesian'
            params.enable_early_stop = true;
            params.enable_targeted_levy = true;
            params.enable_multi_subpop = true;
            params.enable_isolation_penalty = true;
            params.enable_smart_stop = true;
            params.enable_bayesian = false;
            params.enable_adaptive = true;

        case 'no_adaptive'
            params.enable_early_stop = true;
            params.enable_targeted_levy = true;
            params.enable_multi_subpop = true;
            params.enable_isolation_penalty = true;
            params.enable_smart_stop = true;
            params.enable_bayesian = true;
            params.enable_adaptive = false;
    end

    if ~isfield(params, 'enable_early_stop')
        params.enable_early_stop = true;
    end
    if ~isfield(params, 'early_stop_min_iter')
        params.early_stop_min_iter = 50;
    end

    if ~isfield(params, 'enable_multi_subpop')
        params.enable_multi_subpop = true;
    end
    if ~isfield(params, 'enable_bayesian')
        params.enable_bayesian = true;
    end

    if ~params.enable_bayesian
        params.G_weights = [0.35, 0.40, 0.25];
        params.subpop_params.sigma0 = [80, 100, 60];
        params.subpop_params.sigma_min = [5, 8, 3];
        params.subpop_params.w_inertia = [0.7, 0.6, 0.8];
        params.K = 30;
    end

    if params.enable_isolation_penalty
        center_point = [500, 500];
        cover_radius = params.cover_radius;
        effective_priorities = zeros(N_User, 1);
        for i = 1:N_User
            dists_to_others = sqrt(sum((User - User(i,:)).^2, 2));
            neighbor_count = sum(dists_to_others <= cover_radius) - 1;
            dist_to_center = norm(User(i,:) - center_point);
            if neighbor_count == 0
                discount_factor = 0.7 + 0.3 * exp(-dist_to_center / 1000);
            elseif neighbor_count < 3
                discount_factor = 0.9 + 0.1 * exp(-dist_to_center / 1000);
            else
                discount_factor = 1.0;
            end
            effective_priorities(i) = priorities(i) * discount_factor;
        end
        eval_priorities = effective_priorities;
    else
        eval_priorities = priorities;
    end

    subpops = initSubpopulations(N_UAV, User, RRH, eval_priorities, params.subpop_params, Ub, Lb, params.cover_radius, params.D_RU);

    if ~params.enable_multi_subpop
        subpops = {subpops{1}};
        params.G_weights = [1.0];
        params.subpop_params.sigma0 = params.subpop_params.sigma0(1,:);
        params.subpop_params.sigma_min = params.subpop_params.sigma_min(1,:);
        params.subpop_params.w_inertia = params.subpop_params.w_inertia(1);
        params.subpop_params.c = params.subpop_params.c(1);
        params.subpop_params.q = params.subpop_params.q(1);
        params.subpop_params.beta = params.subpop_params.beta(1);
        n_subpops = 1;
    else
        n_subpops = 3;
    end

    if isfield(params, 'K')
        params.K = round(params.K);
        params.K = max(10, min(60, params.K));
    end
    mem_matrix = cell(n_subpops,1);
    for g = 1:n_subpops
        mem_matrix{g} = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, ...
            params.D_UU, params.D_RU);
    end

    cg_curve = zeros(1, params.FES_max);
    curr_curve = zeros(1, params.FES_max);
    energy_consumption = zeros(1, params.FES_max);
    E_remaining_history = zeros(params.FES_max, N_UAV);
    E_remaining = params.E_max * ones(N_UAV, 1);
    weighted_best_curve = zeros(1, params.FES_max);

    pareto_archive = struct('UAV_pos', {}, 'Coverage', {}, 'Energy', {});

    stagnation_counter = zeros(n_subpops, 1);
    prev_fits = zeros(n_subpops, 1);
    last_improve_iter = 1;
    stagnation_generations = 0;

    convergence_window = 10;
    convergence_threshold = 0.2;
    fit_history = zeros(convergence_window, 1);

    decline_window = 5;
    decline_counter = 0;

    capturability_g = zeros(n_subpops, 1);
    for g = 1:n_subpops
        capturability_g(g) = calcCapturability(subpops{g}, 1, params.FES_max, g);
    end
    if ~isfield(params, 'RRH_radius'); params.RRH_radius = 120; end

    [initial_fit, initial_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
        User, eval_priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
        N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
    best_fit = initial_fit;
    bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, eval_priorities, ...
        E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
    if n_subpops >= 2
        [~, weighted_initial] = calcFitness(bestUAV, User, eval_priorities, ...
            E_remaining, params.E_max, params.k_move, 2, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g(2), N_RRH, RRH_type, UAV_type, params);
    else
        weighted_initial = initial_fit;
    end
    cg_curve(1) = best_fit;
    curr_curve(1) = initial_fit;
    weighted_best_curve(1) = weighted_initial;
    energy_consumption(1) = initial_energy;
    E_remaining_history(1,:) = E_remaining';

    initial_cov_high = calcCoverageWithRRH(bestUAV, User(eval_priorities>=3,:), params.cover_radius, RRH, params.RRH_radius);
    initial_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);

    for g = 1:n_subpops
        subpop_fits = zeros(1, size(mem_matrix{g},1));
        for i = 1:size(mem_matrix{g},1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                candidate = reshape(candidate, N_UAV, 2);
            end
            [subpop_fits(i), ~] = calcFitness(candidate, User, eval_priorities, ...
                E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
        end
        prev_fits(g) = max(subpop_fits);
    end

    priorities_sum = sum(eval_priorities);

    % === 去中心化领导 - 各子种群找自己的Leader ===
    local_bests = zeros(3, N_UAV, 2);
    for g_idx = 1:3
        best_f_local = -inf;
        for i_idx = 1:size(mem_matrix{g_idx}, 1)
            cand_temp = squeeze(mem_matrix{g_idx}(i_idx,:,:));
            if size(cand_temp,1)==1; cand_temp = reshape(cand_temp, N_UAV, 2); end
            [f_val, ~] = calcFitness(cand_temp, User, eval_priorities, E_remaining, params.E_max, params.k_move, g_idx, params.subpop_params, N_UAV, params.cover_radius, RRH, capturability_g(g_idx), N_RRH, RRH_type, UAV_type, params);
            if f_val > best_f_local
                best_f_local = f_val;
                local_bests(g_idx, :, :) = cand_temp;
            end
        end
    end

    prev_archive_size = length(pareto_archive);
    mo_stagnation = 0;

    fprintf('初始化完成 (Variant: %s)：适应度=%.4f, 高优先级覆盖率=%.2f%%, 全局覆盖率=%.2f%%\n', ...
        variant, initial_fit, initial_cov_high*100, initial_cov_total*100);

    for iter = 2:params.FES_max
        t = 1 - iter/params.FES_max;

        for g = 1:n_subpops
            capturability_g(g) = calcCapturability(subpops{g}, iter, params.FES_max, g);
        end

        for g = 1:n_subpops
            candidates_init = sampleCandidates(subpops{g}, params.K, N_UAV, Ub, Lb, RRH, ...
                params.D_UU, params.D_RU);
            candidates = zeros(params.K, N_UAV, 2);

            current_mem_size = size(mem_matrix{g}, 1);
            if current_mem_size < params.K
                additional_needed = params.K - current_mem_size;
                additional_candidates = sampleCandidates(subpops{g}, additional_needed, N_UAV, Ub, Lb, RRH, ...
                    params.D_UU, params.D_RU);
                mem_matrix{g} = cat(1, mem_matrix{g}, additional_candidates);
            elseif current_mem_size > params.K
                mem_matrix{g} = mem_matrix{g}(1:params.K, :, :);
            end

            for i = 1:params.K
                cand_i = squeeze(candidates_init(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                    cand_i = cand_i';
                end
                X_mean_g = mean(cand_i, 1);

                for uav_idx = 1:N_UAV
                    X_init = cand_i(uav_idx, :);
                    pos = X_init;

                    mem_idx = min(i, size(mem_matrix{g}, 1));
                    mem_candidate = squeeze(mem_matrix{g}(mem_idx, :, :));
                    if size(mem_candidate, 1) == 2 && size(mem_candidate, 2) == N_UAV
                        mem_candidate = mem_candidate';
                    end
                    mem_ref_pos = mem_candidate(uav_idx, :);

                    if rand >= params.subpop_params.q(g)
                        pos = goaUShape(pos, subpops{g}, mem_ref_pos, t, X_init, g);
                    else
                        pos = goaVShape(pos, subpops{g}, mem_ref_pos, t, X_init, X_mean_g, g);
                    end

                    valid_pos = true;
                    for rrh_idx = 1:N_RRH
                        if norm(pos - RRH(rrh_idx,:)) < params.D_RU
                            valid_pos = false;
                            break;
                        end
                    end
                    if valid_pos
                        for other_uav = 1:N_UAV
                            if other_uav ~= uav_idx
                                other_pos = cand_i(other_uav, :);
                                if norm(pos - other_pos) < params.D_UU
                                    valid_pos = false;
                                    break;
                                end
                            end
                        end
                    end
                    if ~valid_pos
                        pos = X_init;
                    end
                    candidates(i, uav_idx, :) = pos(:)';
                end
            end

            is_stagnant = false;
            if iter > 10 && stagnation_generations > 5
                is_stagnant = true;
            end

            current_mem_size = size(mem_matrix{g}, 1);
            if current_mem_size < params.K
                additional_needed = params.K - current_mem_size;
                additional_candidates = sampleCandidates(subpops{g}, additional_needed, N_UAV, Ub, Lb, RRH, ...
                    params.D_UU, params.D_RU);
                mem_matrix{g} = cat(1, mem_matrix{g}, additional_candidates);
            elseif current_mem_size > params.K
                mem_matrix{g} = mem_matrix{g}(1:params.K, :, :);
            end

            for i = 1:params.K
                cand_i = squeeze(candidates(i, :, :));
                if size(cand_i, 1) == 2 && size(cand_i, 2) == N_UAV
                    cand_i = cand_i';
                end

                for uav_idx = 1:N_UAV
                    pos = cand_i(uav_idx, :);

                    % === 去中心化领导 ===
                    subpop_best_uav = squeeze(local_bests(g, uav_idx, :))';

                    if is_stagnant
                        adjusted_c = params.subpop_params.c(g) * 1.2;
                    else
                        adjusted_c = params.subpop_params.c(g);
                    end

                    % ====== 真实物理空间剩余电量 ======
                    current_fly_dist = norm(pos - [500, 500]);
                    real_E_remaining = params.E_max - params.k_move * current_fly_dist;
                    E_ratio = real_E_remaining / params.E_max;

                    if capturability_g(g) >= adjusted_c
                        pos = goaTurn(pos, subpop_best_uav, capturability_g(g), t);
                        behavior_types{i} = 'turn';
                    else
                        if E_ratio > 0.85
                            if true && is_stagnant
                                high_pri_users = User(eval_priorities >= 3, :);
                                uncovered_targets = [];
                                if ~isempty(high_pri_users)
                                    for ui = 1:size(high_pri_users, 1)
                                        dists_to_uavs = sqrt(sum((bestUAV - high_pri_users(ui, :)).^2, 2));
                                        if min(dists_to_uavs) > params.cover_radius
                                            uncovered_targets = [uncovered_targets; high_pri_users(ui, :)];
                                        end
                                    end
                                end
                                if isempty(uncovered_targets)
                                    if ~isempty(high_pri_users)
                                        uncovered_targets = high_pri_users(randperm(size(high_pri_users, 1), min(3, size(high_pri_users, 1))), :);
                                    else
                                        uncovered_targets = Lb + (Ub - Lb) .* rand(3, 2);
                                    end
                                end
                                % ====== 废除弹弓外推 ======
                                pos = subpop_best_uav + 30 * randn(1, 2);
                                pos = max(Lb, min(Ub, pos));
                            else
                                pos = goaLevy(pos, subpop_best_uav, t);
                            end
                        else
                            pos = goaTurn(pos, subpop_best_uav, capturability_g(g), t);
                        end
                        behavior_types{i} = 'levy';
                    end

                    valid_pos = true;
                    for rrh_idx = 1:N_RRH
                        if norm(pos - RRH(rrh_idx,:)) < params.D_RU
                            valid_pos = false;
                            break;
                        end
                    end
                    if valid_pos
                        for other_uav = 1:N_UAV
                            if other_uav ~= uav_idx
                                other_pos = cand_i(other_uav, :);
                                if norm(pos - other_pos) < params.D_UU
                                    valid_pos = false;
                                    break;
                                end
                            end
                        end
                    end
                    if ~valid_pos
                        pos = cand_i(uav_idx, :);
                    end
                    candidates(i, uav_idx, :) = pos(:)';
                end

                candidate_pos = squeeze(candidates(i, :, :));
                if size(candidate_pos, 1) == 2 && size(candidate_pos, 2) == N_UAV
                    candidate_pos = candidate_pos';
                end
                valid_candidate = true;
                for uav_a = 1:N_UAV
                    for uav_b = uav_a+1:N_UAV
                        if norm(candidate_pos(uav_a,:) - candidate_pos(uav_b,:)) < params.D_UU
                            valid_candidate = false;
                            break;
                        end
                    end
                    if ~valid_candidate
                        break;
                    end
                end
                if ~valid_candidate
                    cand_init_i = squeeze(candidates_init(i, :, :));
                    if size(cand_init_i, 1) == 2 && size(cand_init_i, 2) == N_UAV
                        cand_init_i = cand_init_i';
                    end
                    candidate_pos = cand_init_i;
                end
                candidates(i, :, :) = reshape(candidate_pos, 1, N_UAV, 2);
            end

            mem_matrix{g} = updateMemory(mem_matrix{g}, candidates, User, eval_priorities, ...
                E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
        end

        local_mus = zeros(n_subpops, N_UAV, 2);
        for g = 1:n_subpops
            for uav_idx = 1:N_UAV
                local_mus(g, uav_idx, :) = mean(squeeze(mem_matrix{g}(:, uav_idx, :)), 1);
            end
        end

        for g = 1:n_subpops
            behavior_type = 'turn';
            if capturability_g(g) < params.subpop_params.c(g)
                behavior_type = 'levy';
            end
            subpops{g} = updateSubpopPV(subpops{g}, mem_matrix{g}, squeeze(local_mus(g,:,:)), ...
                params.subpop_params, g, iter, params.FES_max, N_UAV, behavior_type);
        end

        curr_fit_mean = 0;
        total_candidates = 0;
        for g = 1:n_subpops
            for i = 1:size(mem_matrix{g},1)
                candidate = squeeze(mem_matrix{g}(i,:,:));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                [fit_val, ~] = calcFitness(candidate, User, eval_priorities, ...
                    E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                    N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
                curr_fit_mean = curr_fit_mean + fit_val;
                total_candidates = total_candidates + 1;
            end
        end
        curr_fit_mean = curr_fit_mean / total_candidates;

        [curr_fit_best, curr_energy] = calcGlobalFitness(mem_matrix, params.G_weights, ...
            User, eval_priorities, E_remaining, params.E_max, params.k_move, params.subpop_params, ...
            N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);

        curr_curve(iter) = curr_fit_best;

        improvement_threshold = max(1.0, priorities_sum * 1e-4);
        improved_this_iter = false;
        if curr_fit_best > best_fit + improvement_threshold
            best_fit = curr_fit_best;
            bestUAV = calcGlobalBest(mem_matrix, params.G_weights, N_UAV, User, eval_priorities, ...
                E_remaining, params.E_max, params.k_move, params.subpop_params, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
            [~, weighted_best] = calcFitness(bestUAV, User, eval_priorities, ...
                E_remaining, params.E_max, params.k_move, 2, params.subpop_params, ...
                N_UAV, params.cover_radius, RRH, capturability_g(2), N_RRH, RRH_type, UAV_type, params);
            improved_this_iter = true;
        else
            weighted_best = weighted_best_curve(iter-1);
        end
        if ~isscalar(best_fit); best_fit = best_fit(1); end
        cg_curve(iter) = best_fit;
        if ~isscalar(weighted_best); weighted_best = weighted_best(1); end
        weighted_best_curve(iter) = weighted_best;
        if ~isscalar(curr_energy); curr_energy = curr_energy(1); end
        energy_consumption(iter) = curr_energy;

        if iter > 10
            if improved_this_iter
                last_improve_iter = iter;
            end
            stagnation_generations = iter - last_improve_iter;

            threshold = max(8, 15 - floor(iter / 40));

            if stagnation_generations > threshold && params.enable_targeted_levy
                high_pri_users = User(eval_priorities >= 3, :);
                uncovered_targets = [];

                if ~isempty(high_pri_users)
                    for ui = 1:size(high_pri_users, 1)
                        dists_to_uavs = sqrt(sum((bestUAV - high_pri_users(ui, :)).^2, 2));
                        if min(dists_to_uavs) > params.cover_radius
                            uncovered_targets = [uncovered_targets; high_pri_users(ui, :)];
                        end
                    end
                end

                if isempty(uncovered_targets)
                    if ~isempty(high_pri_users)
                        uncovered_targets = high_pri_users(randperm(size(high_pri_users, 1), min(3, size(high_pri_users, 1))), :);
                    else
                        uncovered_targets = Lb + (Ub - Lb) .* rand(3, 2);
                    end
                end

                for g = 1:n_subpops
                    % ================= 核心修复：豁免 G3 =================
                    if g == 3
                        continue;
                    end
                    % ===================================================
                    boost_factor = 1.0 + 0.5 * min(1.0, stagnation_generations / 20);
                    subpops{g}.sigma = min(subpops{g}.sigma * boost_factor, 10 * mean(params.subpop_params.sigma0(:)));

                    if rand < 0.5
                        perturb_count = max(1, round(size(mem_matrix{g}, 1) * 0.15));
                        perturb_idx = randperm(size(mem_matrix{g}, 1), perturb_count);

                        for idx = perturb_idx
                            temp_pos = reshape(mem_matrix{g}(idx, :, :), N_UAV, 2);

                            uavs_to_jump = randperm(N_UAV, randi([1, 2]));
                            for u = uavs_to_jump
                                u_fly_dist = norm(temp_pos(u,:) - [500, 500]);
                                real_E_remaining_u = params.E_max - params.k_move * u_fly_dist;
                                if real_E_remaining_u > 0.85 * params.E_max
                                    target_c = uncovered_targets(randi(size(uncovered_targets, 1)), :);
                                    temp_pos(u, :) = target_c + 15 * randn(1, 2);
                                end
                            end

                            temp_pos = max(Lb, min(Ub, temp_pos));
                            mem_matrix{g}(idx, :, :) = reshape(temp_pos, 1, N_UAV, 2);
                        end
                    end
                end
            end
        end

        curr_subpop_fits = zeros(n_subpops, 1);
        for g = 1:n_subpops
            subpop_fits = zeros(1, size(mem_matrix{g},1));
            for i = 1:size(mem_matrix{g},1)
                candidate = squeeze(mem_matrix{g}(i,:,:));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                [subpop_fits(i), ~] = calcFitness(candidate, User, eval_priorities, ...
                    E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                    N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
            end
            curr_subpop_fits(g) = max(subpop_fits);
        end

        improvement_threshold_subpop = max(1.0, priorities_sum * 1e-4);

        for g = 1:n_subpops
            if curr_subpop_fits(g) > prev_fits(g) + improvement_threshold_subpop
                stagnation_counter(g) = 0;
                prev_fits(g) = curr_subpop_fits(g);
            else
                stagnation_counter(g) = stagnation_counter(g) + 1;
            end

            if iter > 20 && stagnation_counter(g) >= 20
                fprintf('[精英迁移] 迭代 %d: 子种群 G%d 连续 %d 代无改进，触发精英迁移\n', ...
                    iter, g, stagnation_counter(g));
                mem_matrix{g} = migrateElite(mem_matrix, g, Ub, Lb, User, eval_priorities, ...
                    E_remaining, params.E_max, params.k_move, params.subpop_params, ...
                    N_UAV, params.cover_radius, RRH, capturability_g, N_RRH, RRH_type, UAV_type, params);
                stagnation_counter(g) = 0;
                subpop_fits_new = zeros(1, size(mem_matrix{g},1));
                for i = 1:size(mem_matrix{g},1)
                    candidate = squeeze(mem_matrix{g}(i,:,:));
                    if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                        candidate = reshape(candidate, N_UAV, 2);
                    end
                    [subpop_fits_new(i), ~] = calcFitness(candidate, User, eval_priorities, ...
                        E_remaining, params.E_max, params.k_move, g, params.subpop_params, ...
                        N_UAV, params.cover_radius, RRH, capturability_g(g), N_RRH, RRH_type, UAV_type, params);
                end
                curr_subpop_fits(g) = max(subpop_fits_new);
                prev_fits(g) = curr_subpop_fits(g);
            end
        end

        E_remaining = params.E_max * ones(N_UAV, 1);
        E_remaining_history(iter,:) = E_remaining';

        for g = 1:n_subpops
            subpops{g}.prev_mu = subpops{g}.mu;
        end

        curr_cov_high = calcCoverageWithRRH(bestUAV, User(eval_priorities>=3,:), params.cover_radius, RRH, params.RRH_radius);
        curr_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);

        if isfield(params, 'enable_adaptive') && params.enable_adaptive
            algorithm_state = struct();
            algorithm_state.convergence_rate = abs(curr_fit_mean - curr_curve(max(1, iter-1))) / (curr_curve(max(1, iter-1)) + 1e-6);
            algorithm_state.diversity = std(curr_subpop_fits) / (mean(curr_subpop_fits) + 1e-6);
            algorithm_state.stagnation_count = max(stagnation_counter);
            algorithm_state.fitness_improvement = curr_fit_best - best_fit;
            algorithm_state.subpop_fits = curr_subpop_fits;
            algorithm_state.coverage_ratio = curr_cov_total;

            all_positions = [];
            for g = 1:n_subpops
                for i = 1:min(5, size(mem_matrix{g},1))
                    candidate = squeeze(mem_matrix{g}(i,:,:));
                    if size(candidate, 1) == 1
                        candidate = reshape(candidate, N_UAV, 2);
                    end
                    all_positions = [all_positions; candidate];
                end
            end
            if size(all_positions, 1) > 1
                pos_range = max(all_positions) - min(all_positions);
                algorithm_state.search_range = norm(pos_range);
            else
                algorithm_state.search_range = 0;
            end

            if isfield(params, 'enable_adaptive') && params.enable_adaptive
                safe_sigmas = cell(n_subpops, 1);
                for g_sync = 1:n_subpops
                    safe_sigmas{g_sync} = subpops{g_sync}.sigma;
                end
                safe_c = params.subpop_params.c;

                params = adaptiveParams(params, algorithm_state, iter, params.FES_max);

                for g_sync = 1:n_subpops
                    subpops{g_sync}.sigma = min(subpops{g_sync}.sigma, safe_sigmas{g_sync} * 1.05);
                end

                if n_subpops >= 3
                    subpops{3}.sigma = min(subpops{3}.sigma, safe_sigmas{3});
                    params.subpop_params.c(3) = min(params.subpop_params.c(3), safe_c(3));
                end

                params.subpop_params.c = min(params.subpop_params.c, safe_c * 1.1);
            end
        end

        fit_history = circshift(fit_history, -1);
        fit_history(end) = curr_fit_mean;

        fly_distances = sqrt(sum((bestUAV - repmat([500, 500], N_UAV, 1)).^2, 2));
        curr_avg_energy = mean(fly_distances) * params.k_move;
        curr_avg_remaining = params.E_max - curr_avg_energy;
        if mod(iter, 5) == 0 || iter <= 10
            fprintf('迭代 %d/%d (Variant: %s), 当前适应度: %.2f (历史最优: %.2f), 剩余能量: %.2f J (平均能耗: %.2f J), 覆盖率: 高优先级=%.2f%%, 全局=%.2f%%\n', ...
                iter, params.FES_max, variant, curr_fit_mean, best_fit, curr_avg_remaining, curr_avg_energy, ...
                curr_cov_high*100, curr_cov_total*100);
        end

        for g = 1:n_subpops
            for i = 1:size(mem_matrix{g}, 1)
                candidate = squeeze(mem_matrix{g}(i, :, :));
                if size(candidate, 1) == 1 && size(candidate, 2) == N_UAV * 2
                    candidate = reshape(candidate, N_UAV, 2);
                end
                cand_cov = calcCoverageWithRRH(candidate, User, params.cover_radius, RRH, params.RRH_radius);
                cov_count = cand_cov * N_User;
                cand_fly_dist = sqrt(sum((candidate - repmat([500, 500], N_UAV, 1)).^2, 2));
                cand_energy = sum(cand_fly_dist * params.k_move);
                pareto_archive = updateParetoArchive(pareto_archive, candidate, cov_count, cand_energy);
            end
        end

        current_archive_size = length(pareto_archive);
        if current_archive_size > prev_archive_size
            mo_stagnation = 0;
        else
            mo_stagnation = mo_stagnation + 1;
        end
        prev_archive_size = current_archive_size;

        if params.enable_smart_stop && params.enable_early_stop && iter > 100
            if mo_stagnation >= 25
                fprintf('\n[提前停止] 迭代 %d/%d: 多目标彻底收敛，连续 %d 代无新 Pareto 解\n', iter, params.FES_max, mo_stagnation);
                break;
            end
        end
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_distances = sqrt(sum((bestUAV - center_point).^2, 2));
    fly_energy_costs = params.k_move * fly_distances;
    final_E_remaining = max(0, params.E_max - fly_energy_costs');

    final_cov_high = calcCoverageWithRRH(bestUAV, User(eval_priorities>=3,:), params.cover_radius, RRH, params.RRH_radius);
    final_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);

    cov_count = final_cov_total * N_User;
    curr_total_energy = sum(fly_energy_costs);
    pareto_archive = updateParetoArchive(pareto_archive, bestUAV, cov_count, curr_total_energy);

    fprintf('\n========== 算法完成 (Variant: %s) ==========\n', variant);
    if iter < params.FES_max
        fprintf('提前停止于迭代 %d/%d\n', iter, params.FES_max);
    else
        fprintf('达到最大迭代次数 %d\n', params.FES_max);
    end

    center_point = [500, 500];
    total_fly_distance = sum(sqrt(sum((bestUAV - repmat(center_point, N_UAV, 1)).^2, 2)));
    total_energy_consumed = params.k_move * total_fly_distance;
    avg_energy_consumed = total_energy_consumed / N_UAV;

    fprintf('最终结果：\n');
    fprintf('  1. 历史最优适应度：%.2f（优先级和）\n', best_fit);
    fprintf('  2. 最终当前适应度：%.2f（优先级和）\n', curr_curve(iter));
    fprintf('  3. 总飞行能耗：%.2f J，平均每UAV：%.2f J\n', total_energy_consumed, avg_energy_consumed);
    fprintf('  4. 高优先级覆盖率：%.2f%%\n', final_cov_high*100);
    fprintf('  5. 全局覆盖率：%.2f%%\n', final_cov_total*100);
    fprintf('  6. 平均每UAV飞行距离：%.2f m\n', total_fly_distance / N_UAV);
    fprintf('============================\n');
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos,1)
        dists_uav = sqrt(sum((UAV_pos - User_pos(i,:)).^2, 2));
        covered_by_uav = any(dists_uav <= UAV_radius);

        if size(RRH,1) > 0
            dists_rrh = sqrt(sum((RRH - User_pos(i,:)).^2, 2));
            covered_by_rrh = any(dists_rrh <= RRH_radius);
        else
            covered_by_rrh = false;
        end

        if covered_by_uav || covered_by_rrh
            covered = covered + 1;
        end
    end
    cov_ratio = covered / size(User_pos,1);
end

function new_pos = goaUShape(pos, subpop, mem_ref_pos, t, X_init, g)
    r2 = 2*pi*rand;
    a_coeffs = [0.6, 0.7, 0.5];
    r4 = rand;
    A_g = (2*r4 - 1) * a_coeffs(g);

    X_r = mem_ref_pos(:)';
    X_init = X_init(:)';

    sigma_mean = mean(subpop.sigma(:));
    new_pos = X_init + 3*cos(r2)*t*sigma_mean + A_g*(X_r - X_init);
end

function new_pos = goaVShape(pos, subpop, mem_ref_pos, t, X_init, X_mean, g)
    r3 = rand;
    x = 2*pi*r3;

    if x >= 0 && x < pi
        V_x = -x/pi + 1;
    else
        V_x = x/pi - 1;
    end

    b_coeffs = [0.5, 0.6, 0.4];
    r5 = rand;
    B_g = (2*r5 - 1) * b_coeffs(g);

    sigma_mean = mean(subpop.sigma(:));
    X_init = X_init(:)';
    X_mean = X_mean(:)';
    new_pos = X_init + 3*V_x*t*sigma_mean + B_g*(X_mean - X_init);
end

function new_pos = goaTurn(pos, global_best_uav, cap, t)
    delta = cap * norm(pos - global_best_uav);
    direction = sign(global_best_uav - pos);
    theta = randn * 0.2;
    rot_matrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    direction = direction * rot_matrix;
    new_pos = pos(:)' + t*delta.*direction;
end

function new_pos = goaLevy(pos, global_best_uav, t)
    beta_levy = 1.5;
    sigma_levy = (gamma(1+beta_levy)*sin(pi*beta_levy/2) / ...
        (gamma((1+beta_levy)/2)*beta_levy*2^((beta_levy-1)/2)))^(1/beta_levy);

    mu = randn(1,2);
    v = randn(1,2);
    levy = 0.05 * (mu * sigma_levy) ./ (abs(v).^(1/beta_levy) + 1e-10);

    new_pos = pos(:)' + t * levy .* sign(global_best_uav(:)' - pos(:)');
    new_pos = max(-1, min(1, new_pos));
end

function archive = updateParetoArchive(archive, new_pos, new_cov, new_energy)
    is_dominated = false;
    indices_to_remove = [];

    for i = 1:length(archive)
        if (archive(i).Coverage >= new_cov && archive(i).Energy <= new_energy) && ...
           (archive(i).Coverage > new_cov || archive(i).Energy < new_energy)
            is_dominated = true;
            break;
        end

        if (new_cov >= archive(i).Coverage && new_energy <= archive(i).Energy) && ...
           (new_cov > archive(i).Coverage || new_energy < archive(i).Energy)
            indices_to_remove = [indices_to_remove, i];
        end
    end

    if ~is_dominated
        archive(indices_to_remove) = [];

        new_entry.UAV_pos = new_pos;
        new_entry.Coverage = new_cov;
        new_entry.Energy = new_energy;
        archive(end+1) = new_entry;
    end
end
