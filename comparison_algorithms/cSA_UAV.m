function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = cSA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
    end

    lamda = 10;
    Np = 300;

    mu = zeros(N_UAV, 2);
    sicma = lamda * ones(N_UAV, 2);

    population = zeros(N_UAV, 2);
    for j = 1:N_UAV
        population(j, 1) = Lb(1) + rand() * (Ub(1) - Lb(1));
        population(j, 2) = Lb(2) + rand() * (Ub(2) - Lb(2));
    end

    if ~checkConstraints(population, params.D_UU, params.D_RU, RRH)
        population = enforceConstraints(population, params.D_UU, params.D_RU, RRH);
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((population - center_point).^2, 2));
    fly_energy = params.k_move * fly_dist;
    E_curr = max(0, params.E_max - fly_energy);

    [best_fit, ~] = calcFitness(population, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
        params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);
    bestUAV = population;

    cg_curve = zeros(1, max_iter);
    cg_curve(1) = best_fit;

    pareto_archive = struct('UAV_pos', {}, 'Coverage', {}, 'Energy', {});

    iter_count = 0;

    while iter_count < max_iter
        a = 2 - iter_count * (2 / max_iter);

        new_population = zeros(N_UAV, 2);
        for j = 1:N_UAV
            new_population(j, 1) = Lb(1) + rand() * (Ub(1) - Lb(1));
            new_population(j, 2) = Lb(2) + rand() * (Ub(2) - Lb(2));
        end

        for i = 1:N_UAV
            for j = 1:2
                r1 = 2 * a * rand() - a;
                r2 = 2 * pi * rand();
                r3 = 2 * rand();
                new_population(i, j) = new_population(i, j) + ...
                    (r1 * sin(r2) * (r3 * bestUAV(i, j) - new_population(i, j)));
            end

            flagub = new_population(i, :) > Ub;
            new_population(i, flagub) = 2 * Ub(flagub) - new_population(i, flagub);
            flaglb = new_population(i, :) < Lb;
            new_population(i, flaglb) = 2 * Lb(flaglb) - new_population(i, flaglb);
        end

        for i = 1:N_UAV
            tmp_UAV = bestUAV;
            tmp_UAV(i, :) = new_population(i, :);
            uavindex = 1:N_UAV;
            uavindex(i) = [];

            if all(sqrt(sum((tmp_UAV(i, :) - tmp_UAV(uavindex, :)).^2, 2)) >= params.D_UU) && ...
               all(sqrt(sum((tmp_UAV(i, :) - RRH(1:end, :)).^2, 2)) >= params.D_RU)

                iter_count = iter_count + 1;

                fly_dist = sqrt(sum((tmp_UAV - center_point).^2, 2));
                fly_energy = params.k_move * fly_dist;
                E_curr = max(0, params.E_max - fly_energy);

                [tmp_fitness, ~] = calcFitness(tmp_UAV, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                    params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);

                winner = 2 * (bestUAV(i, :) - Lb) ./ (Ub - Lb) - 1;
                loser = 2 * (new_population(i, :) - Lb) ./ (Ub - Lb) - 1;

                if tmp_fitness > best_fit
                    winner = 2 * (new_population(i, :) - Lb) ./ (Ub - Lb) - 1;
                    loser = 2 * (bestUAV(i, :) - Lb) ./ (Ub - Lb) - 1;
                    bestUAV(i, :) = new_population(i, :);
                    best_fit = tmp_fitness;
                end

                for k = 1:2
                    mut = mu(i, k);
                    mu(i, k) = mut + (1 / Np) * (winner(k) - loser(k));
                    tt = sicma(i, k)^2 + mut^2 - mu(i, k)^2 + (1 / Np) * (winner(k)^2 - loser(k)^2);
                    if tt > 0
                        sicma(i, k) = sqrt(tt);
                    else
                        sicma(i, k) = 10;
                    end
                end

                if iter_count > max_iter
                    break;
                end

                cg_curve(iter_count) = best_fit;

                cov = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);
                pareto_archive = updateParetoArchive(pareto_archive, bestUAV, cov * N_User, sum(fly_energy));
            end
        end

        if mod(iter_count, 50) == 0
            fprintf('cSA iter %d/%d, Best fitness: %.4f\n', iter_count, max_iter, best_fit);
        end
    end

    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
    best_energy = sum(params.k_move * fly_dist);
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
