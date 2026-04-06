function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = PSO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    pop_size = 40;
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
        pop_size = 40;
    end

    w = 0.4;
    c1 = 1.5;
    c2 = 1.5;

    n_vars = N_UAV * 2;
    initial_pos = ones(N_UAV, 1) * [500, 500];

    population = zeros(pop_size, n_vars);
    velocities = zeros(pop_size, n_vars);
    pbest = zeros(pop_size, n_vars);
    pbest_obj = zeros(pop_size, 1);
    pbest_energy = zeros(pop_size, N_UAV);

    E_remaining = params.E_max * ones(N_UAV, 1);

    center_point = [500, 500];
    jitter = 10;

    for i = 1:pop_size
        for j = 1:N_UAV
            init_x = center_point(1) + jitter * randn();
            init_y = center_point(2) + jitter * randn();
            population(i, (j-1)*2+1) = max(Lb(1), min(Ub(1), init_x));
            population(i, (j-1)*2+2) = max(Lb(2), min(Ub(2), init_y));
        end
        velocities(i, :) = (rand(1, n_vars) - 0.5) * (Ub(1) - Lb(1)) * 0.1;
    end

    for i = 1:pop_size
        uav_pos = reshape(population(i, :), N_UAV, 2);
        if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
            uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
            population(i, :) = reshape(uav_pos, 1, N_UAV * 2);
        end
        center_point = repmat([500, 500], N_UAV, 1);
        fly_dist = sqrt(sum((uav_pos - center_point).^2, 2));
        fly_energy = params.k_move * fly_dist;
        E_curr = max(0, params.E_max - fly_energy);
        [fitness, ~] = calcFitness(uav_pos, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
            params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
        pbest(i, :) = population(i, :);
        pbest_obj(i) = fitness;
        pbest_energy(i, :) = E_curr';
    end

    [best_fit, best_idx] = max(pbest_obj);
    bestUAV = reshape(pbest(best_idx, :), N_UAV, 2);
    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
    best_energy = sum(params.k_move * fly_dist);

    cg_curve = zeros(1, max_iter);
    cg_curve(1) = best_fit;

    pareto_archive = struct('UAV_pos', {}, 'Coverage', {}, 'Energy', {});

    for iter = 2:max_iter
        gbest = bestUAV(:)';

        for i = 1:pop_size
            r1 = rand(1, n_vars);
            r2 = rand(1, n_vars);
            velocities(i, :) = w * velocities(i, :) + ...
                c1 * r1 .* (pbest(i, :) - population(i, :)) + ...
                c2 * r2 .* (gbest - population(i, :));

            population(i, :) = population(i, :) + velocities(i, :);

            dim_idx = 1;
            for j = 1:N_UAV
                population(i, dim_idx) = max(Lb(1), min(Ub(1), population(i, dim_idx)));
                population(i, dim_idx+1) = max(Lb(2), min(Ub(2), population(i, dim_idx+1)));
                dim_idx = dim_idx + 2;
            end

            uav_pos = reshape(population(i, :), N_UAV, 2);
            if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
                uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
                population(i, :) = reshape(uav_pos, 1, N_UAV * 2);
            end

            center_point = repmat([500, 500], N_UAV, 1);
            fly_dist = sqrt(sum((uav_pos - center_point).^2, 2));
            fly_energy = params.k_move * fly_dist;
            E_curr = max(0, params.E_max - fly_energy);

            [fitness, ~] = calcFitness(uav_pos, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);

            if fitness > pbest_obj(i)
                pbest(i, :) = population(i, :);
                pbest_obj(i) = fitness;
                pbest_energy(i, :) = E_curr';
            end

            if fitness > best_fit
                best_fit = fitness;
                bestUAV = uav_pos;
                fly_dist = sqrt(sum((uav_pos - center_point).^2, 2));
                best_energy = sum(params.k_move * fly_dist);
            end

            cov = calcCoverageWithRRH(uav_pos, User, params.cover_radius, RRH, params.RRH_radius);
            pareto_archive = updateParetoArchive(pareto_archive, uav_pos, cov * N_User, sum(fly_energy));
        end

        cg_curve(iter) = best_fit;

        if mod(iter, 50) == 0
            fprintf('PSO iter %d/%d, Best fitness: %.4f\n', iter, max_iter, best_fit);
        end
    end
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

function feasible = checkConstraints(uav_pos, D_UU, D_RU, RRH)
    N_UAV = size(uav_pos, 1);
    feasible = true;

    for i = 1:N_UAV
        for j = i+1:N_UAV
            dist = sqrt(sum((uav_pos(i, :) - uav_pos(j, :)).^2, 2));
            if dist < D_UU
                feasible = false;
                return;
            end
        end
    end

    for i = 1:N_UAV
        for j = 1:size(RRH, 1)
            dist = sqrt(sum((uav_pos(i, :) - RRH(j, :)).^2, 2));
            if dist < D_RU
                feasible = false;
                return;
            end
        end
    end
end

function uav_pos = enforceConstraints(uav_pos, D_UU, D_RU, RRH)
    N_UAV = size(uav_pos, 1);
    max_iter = 100;

    for iter = 1:max_iter
        changed = false;

        for i = 1:N_UAV
            for j = i+1:N_UAV
                dist = sqrt(sum((uav_pos(i, :) - uav_pos(j, :)).^2, 2));
                if dist < D_UU && dist > 0
                    direction = (uav_pos(i, :) - uav_pos(j, :)) / dist;
                    uav_pos(i, :) = uav_pos(i, :) + direction * (D_UU - dist) * 0.5;
                    uav_pos(j, :) = uav_pos(j, :) - direction * (D_UU - dist) * 0.5;
                    changed = true;
                end
            end
        end

        for i = 1:N_UAV
            for j = 1:size(RRH, 1)
                dist = sqrt(sum((uav_pos(i, :) - RRH(j, :)).^2, 2));
                if dist < D_RU && dist > 0
                    direction = (uav_pos(i, :) - RRH(j, :)) / dist;
                    uav_pos(i, :) = uav_pos(i, :) + direction * (D_RU - dist);
                    changed = true;
                end
            end
        end

        if ~changed
            break;
        end
    end
end
