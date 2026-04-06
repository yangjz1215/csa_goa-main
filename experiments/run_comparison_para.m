function results = run_comparison_para(varargin)
    fprintf('========== 对比实验开始 (并行版本) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'comparison_algorithms')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 30);
    addParameter(p, 'map_name', 'Map1_Medium');
    addParameter(p, 'verbose', false);
    addParameter(p, 'n_workers', 3);
    parse(p, varargin{:});
    n_runs = p.Results.n_runs;
    map_name = p.Results.map_name;
    n_workers = p.Results.n_workers;

    maps_dir = fullfile(project_dir, 'maps');
    map_file = fullfile(maps_dir, [map_name, '.mat']);

    if ~exist(map_file, 'file')
        error('地图文件不存在: %s', map_file);
    end

    fprintf('--- 加载固定地图: %s ---\n', map_name);
    map_data = load(map_file);
    User = map_data.User;
    priorities = map_data.priorities;
    N_User = map_data.N_User;

    N_User_actual = size(User, 1);
    if N_User_actual == 200
        N_UAV = 8;
        config_suffix = 'Small';
    elseif N_User_actual == 500
        N_UAV = 15;
        config_suffix = 'Medium';
    else
        N_UAV = 25;
        config_suffix = 'Large';
    end

    fprintf('  用户数: %d, UAV数: %d\n', N_User, N_UAV);

    config_file = fullfile(maps_dir, ['Map_', config_suffix, '_Config.mat']);
    if exist(config_file, 'file')
        config_data = load(config_file);
        RRH = config_data.RRH;
        RRH_type = config_data.RRH_type;
        N_RRH = config_data.N_RRH;
        N_eRRH = config_data.N_eRRH;
        UAV_type = config_data.UAV_type;
        params.D = config_data.D;
        params.C = config_data.C;
        params.DT = config_data.DT;
    else
        Lb = [0, 0];
        Ub = [1000, 1000];
        N_RRH = 10;
        N_eRRH = 4;
        RRH = Lb + (Ub - Lb) .* rand(N_RRH, 2);
        RRH_type = zeros(N_RRH, 1);
        RRH_type(1:N_eRRH) = 1;
        UAV_type = zeros(N_UAV, 1);
        UAV_type(1:floor(N_UAV * 0.3)) = 1;
        params.D = [];
        params.C = [];
        params.DT = [];
    end

    Lb = [0, 0];
    Ub = [1000, 1000];

    params.D_UU = 10;
    params.D_RU = 10;
    params.cover_radius = 150;
    params.RRH_radius = 150;
    params.E_max = 50000;
    params.k_move = 15;
    params.bandwidth = 1e6;

    params.Pho = 100;
    params.ki = 1e-27;
    params.PtxU = 10;
    params.PtxEU = 5;
    params.Ptx = 1;
    params.PtxR = 10;
    params.alpha0 = 1.42e-4;
    params.sigma2 = 3.98e-12;
    params.f_eUAV = 2e9;
    params.f_eRRH = 2e9;
    params.f_BBU = 4e9;

    params.FES_max = 300;
    params.K = 40;

    params.G_weights = [0.4, 0.3, 0.3];
    params.subpop_params = struct();
    params.subpop_params.mu0 = [500, 500];
    params.subpop_params.sigma0 = [150, 150; 120, 120; 80, 80];
    params.subpop_params.sigma_min = [5, 8, 3];
    params.subpop_params.w_inertia = [0.7, 0.6, 0.8];
    params.subpop_params.c = [0.15, 0.10, 0.08];
    params.subpop_params.q = [0.6, 0.5, 0.4];
    params.subpop_params.beta = [0.8, 0.7, 0.6];

    algorithms = {
        'cSA_GOA', 'cSA-GOA (Proposed)';
        'PSO', 'PSO (Particle Swarm Optimization)';
        'GA', 'GA (Genetic Algorithm)';
        'GOA', 'GOA (Grasshopper Optimization)';
        'cSA', 'cSA (Compact Sine Algorithm)';
        'GWO', 'GWO (Grey Wolf Optimizer)'
    };

    reference_point = [1.0, 100000];

    if isempty(gcp('nocreate'))
        fprintf('启动并行池 (%d workers)...\n', n_workers);
        parpool('local', n_workers);
    end

    results = struct();
    results.map_name = map_name;
    results.N_User = N_User;
    results.N_UAV = N_UAV;

    for alg_idx = 1:size(algorithms, 1)
        alg_name = algorithms{alg_idx, 1};
        alg_desc = algorithms{alg_idx, 2};
        fprintf('\n--- 运行算法: %s (并行) ---\n', alg_desc);

        best_fits = zeros(n_runs, 1);
        energies = zeros(n_runs, 1);
        cov_high = zeros(n_runs, 1);
        cov_total = zeros(n_runs, 1);
        iter_counts = zeros(n_runs, 1);
        convergence_curves = cell(n_runs, 1);
        hv_values = zeros(n_runs, 1);
        igd_values = zeros(n_runs, 1);
        spread_values = zeros(n_runs, 1);
        pareto_sizes = zeros(n_runs, 1);
        pareto_fronts_cell = cell(n_runs, 1);

        parfor run = 1:n_runs
            stream = RandStream('mt19937ar', 'Seed', run * 200 + alg_idx * 1000);
            RandStream.setGlobalStream(stream);

            best_fit = 0;
            bestUAV = zeros(N_UAV, 2);
            cg_curve = zeros(1, 300);
            pareto_archive = struct('Coverage', {}, 'Energy', {}, 'UAV_pos', {});

            switch alg_name
                case 'cSA_GOA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        model0_proposed(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'PSO'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        PSO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        GA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GOA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        GOA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'cSA'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        cSA_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'GWO'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        GWO_UAV(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
            end

            if ~isempty(pareto_archive) && length(pareto_archive) > 1
                arch_cov = [pareto_archive.Coverage];
                arch_energy = [pareto_archive.Energy];

                norm_cov = (arch_cov - min(arch_cov)) / (max(arch_cov) - min(arch_cov) + 1e-6);
                norm_eng = (arch_energy - min(arch_energy)) / (max(arch_energy) - min(arch_energy) + 1e-6);
                distances_to_ideal = sqrt((1 - norm_cov).^2 + (0 - norm_eng).^2);
                [~, idx_knee] = min(distances_to_ideal);

                best_fits(run) = best_fit;
                energies(run) = pareto_archive(idx_knee).Energy;

                if isfield(pareto_archive(idx_knee), 'UAV_pos')
                    bestUAV = pareto_archive(idx_knee).UAV_pos;
                end
            else
                best_fits(run) = best_fit;
                center_point = repmat([500, 500], N_UAV, 1);
                fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
                energies(run) = sum(params.k_move * fly_dist);
            end

            cov_high(run) = calcCoverageWithRRH(bestUAV, User(priorities>=3,:), params.cover_radius, RRH, params.RRH_radius) * 100;
            cov_total(run) = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius) * 100;

            iter_counts(run) = length(cg_curve);
            convergence_curves{run} = cg_curve;

            if ~isempty(pareto_archive) && length(pareto_archive) > 1
                pareto_front = zeros(length(pareto_archive), 2);
                for p_idx = 1:length(pareto_archive)
                    pareto_front(p_idx, 1) = pareto_archive(p_idx).Coverage / N_User;
                    pareto_front(p_idx, 2) = pareto_archive(p_idx).Energy;
                end
                pareto_sizes(run) = length(pareto_archive);
                pareto_fronts_cell{run} = pareto_front;

                try
                    % HV量纲归一化：能耗目标归一化到[0,1]区间
                    % 参考点是[1.0, 100000]，覆盖率已是0~1，能耗需除以100000
                    ref_point_norm = [1.0, 1.0];
                    pareto_front_norm = pareto_front;
                    pareto_front_norm(:, 2) = pareto_front(:, 2) / 100000;

                    metrics = calculate_all_metrics(pareto_front_norm, [], ref_point_norm);
                    hv_values(run) = metrics.hv;
                    spread_values(run) = metrics.spread;
                catch
                    hv_values(run) = 0;
                    spread_values(run) = NaN;
                end
            else
                pareto_sizes(run) = 0;
                hv_values(run) = 0;
                spread_values(run) = NaN;
                pareto_fronts_cell{run} = [];
            end
        end

        pareto_fronts{alg_idx} = pareto_fronts_cell;

        results.(alg_name) = struct();
        results.(alg_name).description = alg_desc;
        results.(alg_name).best_fits = best_fits;
        results.(alg_name).energies = energies;
        results.(alg_name).cov_high = cov_high;
        results.(alg_name).cov_total = cov_total;
        results.(alg_name).iter_counts = iter_counts;
        results.(alg_name).convergence_curves = convergence_curves;
        results.(alg_name).mean_fitness = mean(best_fits);
        results.(alg_name).std_fitness = std(best_fits);
        results.(alg_name).mean_energy = mean(energies);
        results.(alg_name).mean_cov_high = mean(cov_high);
        results.(alg_name).mean_cov_total = mean(cov_total);
        results.(alg_name).hv_values = hv_values;
        results.(alg_name).mean_hv = mean(hv_values);
        results.(alg_name).std_hv = std(hv_values);
        results.(alg_name).igd_values = igd_values;
        results.(alg_name).mean_igd = mean(igd_values);
        results.(alg_name).spread_values = spread_values;
        results.(alg_name).mean_spread = mean(spread_values);
        results.(alg_name).mean_pareto_size = mean(pareto_sizes);

        fprintf('  >> %s 平均结果:\n', alg_desc);
        fprintf('     Mean Fitness: %.2f +/- %.2f\n', mean(best_fits), std(best_fits));
        fprintf('     Mean Energy: %.2f J\n', mean(energies));
        fprintf('     Mean High-Priority Coverage: %.2f%%\n', mean(cov_high));
        fprintf('     Mean HV: %.4f +/- %.4f\n', mean(hv_values), std(hv_values));
    end

    fprintf('\n========== 计算IGD (使用合并Pareto前沿作为参考) ==========\n');
    all_pareto_points = [];
    for alg_idx = 1:size(algorithms, 1)
        for run = 1:n_runs
            pf = pareto_fronts{alg_idx}{run};
            if ~isempty(pf)
                all_pareto_points = [all_pareto_points; pf];
            end
        end
    end

    if ~isempty(all_pareto_points)
        true_front = extractNonDominated(all_pareto_points);
        fprintf('  合并Pareto解数量: %d, 非支配解数量: %d\n', size(all_pareto_points, 1), size(true_front, 1));

        % 归一化true_front（用于IGD计算）
        true_front_norm = true_front;
        true_front_norm(:, 2) = true_front(:, 2) / 100000;

        for alg_idx = 1:size(algorithms, 1)
            alg_name = algorithms{alg_idx, 1};
            igd_values = zeros(n_runs, 1);
            for run = 1:n_runs
                pf = pareto_fronts{alg_idx}{run};
                if ~isempty(pf)
                    % 归一化前沿用于IGD计算
                    pf_norm = pf;
                    pf_norm(:, 2) = pf(:, 2) / 100000;
                    igd_values(run) = igd(pf_norm, true_front_norm);
                else
                    igd_values(run) = NaN;
                end
            end
            results.(alg_name).igd_values = igd_values;
            results.(alg_name).mean_igd = mean(igd_values);
            results.(alg_name).std_igd = std(igd_values);
        end
    end

    fprintf('\n========== 对比实验完成 (地图: %s) ==========\n', map_name);
    fprintf('\n========== 结果汇总表格 (多目标三剑客指标) ==========\n');
    fprintf('%-18s | %-8s | %-8s | %-8s | %-8s | %-10s | %-10s | %-8s\n', ...
        'Algorithm', 'Fitness', 'Energy(J)', 'HighPri%', 'Total%', 'HV', 'IGD', 'Spread');
    fprintf('%s\n', repmat('-', 1, 115));
    for alg_idx = 1:size(algorithms, 1)
        alg_name = algorithms{alg_idx, 1};
        r = results.(alg_name);
        fprintf('%-18s | %-8.2f | %-8.2f | %-8.2f | %-8.2f | %-8.4f±%-6.4f | %-8.4f±%-6.4f | %-8.4f\n', ...
            alg_name, r.mean_fitness, r.mean_energy, ...
            r.mean_cov_high, r.mean_cov_total, ...
            r.mean_hv, r.std_hv, r.mean_igd, r.std_igd, r.mean_spread);
    end
    fprintf('================================\n');
    fprintf('注: HV↑越大越好 | IGD↓越小越好 | Spread↑分布越均匀\n');

    results_file = fullfile(project_dir, 'experiments', ...
        ['comparison_results_para_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results');
    fprintf('结果已保存: %s\n', results_file);
end

function cov_ratio = calcCoverageWithRRH(UAV_pos, User_pos, UAV_radius, RRH, RRH_radius)
    covered = 0;
    for i = 1:size(User_pos, 1)
        user = User_pos(i, :);
        dist_UAV = min(sqrt(sum((UAV_pos - repmat(user, size(UAV_pos, 1), 1)).^2, 2)));
        dist_RRH = min(sqrt(sum((RRH - repmat(user, size(RRH, 1), 1)).^2, 2)));
        if dist_UAV <= UAV_radius || dist_RRH <= RRH_radius
            covered = covered + 1;
        end
    end
    cov_ratio = covered / size(User_pos, 1);
end

function pf = extractNonDominated(points)
    n = size(points, 1);
    is_dominated = false(n, 1);
    for i = 1:n
        for j = 1:n
            if i ~= j
                if points(j, 1) >= points(i, 1) && points(j, 2) <= points(i, 2)
                    if points(j, 1) > points(i, 1) || points(j, 2) < points(i, 2)
                        is_dominated(i) = true;
                        break;
                    end
                end
            end
        end
    end
    pf = points(~is_dominated, :);
end