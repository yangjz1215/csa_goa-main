function [best_fit, bestUAV, cg_curve, best_energy, pareto_archive] = GOA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities)
    pop_size = 40;
    max_iter = 300;
    if isfield(params, 'FES_max')
        max_iter = params.FES_max;
        pop_size = 40;
    end

    n_vars = N_UAV * 2;
    dim = n_vars;

    population = zeros(pop_size, n_vars);

    center_point = [500, 500];
    jitter = 10;

    for i = 1:pop_size
        for j = 1:N_UAV
            init_x = center_point(1) + jitter * randn();
            init_y = center_point(2) + jitter * randn();
            population(i, (j-1)*2+1) = max(Lb(1), min(Ub(1), init_x));
            population(i, (j-1)*2+2) = max(Lb(2), min(Ub(2), init_y));
        end
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
        fitness_pop(i) = fitness;
    end

    [best_fit, best_idx] = max(fitness_pop);
    bestUAV = reshape(population(best_idx, :), N_UAV, 2);
    center_point = repmat([500, 500], N_UAV, 1);
    fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
    best_energy = sum(params.k_move * fly_dist);

    cg_curve = zeros(1, max_iter);
    cg_curve(1) = best_fit;

    pareto_archive = struct('UAV_pos', {}, 'Coverage', {}, 'Energy', {});

    for iter = 1:max_iter
        t = 1 - iter / max_iter;
        t2 = 1 + iter / max_iter;
        a = 2 * cos(2 * pi * rand) * t;
        b = 2 * Vshape(2 * pi * rand) * t;
        A = (2 * rand - 1) * a;
        B = (2 * rand - 1) * b;

        if rand > 0.5
            for i = 1:pop_size
                if rand >= 0.5
                    u1 = unifrnd(-a, a, 1, dim);
                    idx = randi(pop_size);
                    u2 = A * (population(i, :) - population(idx, :));
                    MX(i, :) = population(i, :) + u1 + u2;
                else
                    v1 = unifrnd(-b, b, 1, dim);
                    v2 = B * (population(i, :) - mean(population, 1));
                    MX(i, :) = population(i, :) + v1 + v2;
                end
            end
        else
            vel = 1.5;
            M = 2.5;
            L = 0.2 + (2 - 0.2) * rand();
            r = M * vel^2 / L;
            Captureability = 1 / (r * t2);

            for i = 1:pop_size
                if Captureability > 0.15
                    delta = Captureability .* abs(population(i, :) - bestUAV(:)');
                    MX(i, :) = t * delta .* (population(i, :) - bestUAV(:)') + population(i, :);
                else
                    P = Levy(dim);
                    MX(i, :) = bestUAV(:)' - (bestUAV(:)' - population(i, :)) .* P * t;
                end
            end
        end

        for i = 1:pop_size
            MX(i, :) = max(Lb(1), min(Ub(1), MX(i, :)));
            MX(i, 2:2:end) = max(Lb(2), min(Ub(2), MX(i, 2:2:end)));

            uav_pos = reshape(MX(i, :), N_UAV, 2);
            if ~checkConstraints(uav_pos, params.D_UU, params.D_RU, RRH)
                uav_pos = enforceConstraints(uav_pos, params.D_UU, params.D_RU, RRH);
                MX(i, :) = reshape(uav_pos, 1, N_UAV * 2);
            end

            center_point = repmat([500, 500], N_UAV, 1);
            fly_dist = sqrt(sum((uav_pos - center_point).^2, 2));
            fly_energy = params.k_move * fly_dist;
            E_curr = max(0, params.E_max - fly_energy);

            [fitness, ~] = calcFitness(uav_pos, User, priorities, E_curr, params.E_max, params.k_move, 1, ...
                params.subpop_params, N_UAV, params.cover_radius, RRH, 0.5, N_RRH, RRH_type, UAV_type, params);

            if fitness > fitness_pop(i)
                fitness_pop(i) = fitness;
                population(i, :) = MX(i, :);
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
            fprintf('GOA iter %d/%d, Best fitness: %.4f\n', iter, max_iter, best_fit);
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

function o = Vshape(v)
    o = abs(v);
end

function step = Levy(d)
    beta = 1.5;
    sigma = (gamma(1 + beta) * sin(pi * beta / 2) / (gamma((1 + beta) / 2) * beta * 2 ^ ((beta - 1) / 2))) ^ (1 / beta);
    u = randn(1, d) * sigma;
    v = randn(1, d);
    step = u ./ abs(v) .^ (1 / beta);
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
