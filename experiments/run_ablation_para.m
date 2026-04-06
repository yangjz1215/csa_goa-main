function results = run_ablation_para(varargin)
    fprintf('========== 消融实验开始 (并行版本) ==========\n');

    script_dir = fileparts(mfilename('fullpath'));
    project_dir = fileparts(script_dir);
    addpath(genpath(project_dir));
    addpath(genpath(fullfile(project_dir, 'ablation')));
    addpath(genpath(fullfile(project_dir, 'performance_metrics')));

    p = inputParser;
    addParameter(p, 'n_runs', 30);
    addParameter(p, 'map_name', 'Map1_Medium');
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
    params.K = 15;

    params.G_weights = [0.4, 0.3, 0.3];
    params.subpop_params = struct();
    params.subpop_params.mu0 = [500, 500];
    params.subpop_params.sigma0 = [150, 150; 120, 120; 80, 80];
    params.subpop_params.sigma_min = [5, 8, 3];
    params.subpop_params.w_inertia = [0.7, 0.6, 0.8];
    params.subpop_params.c = [0.15, 0.10, 0.08];
    params.subpop_params.q = [0.6, 0.5, 0.4];
    params.subpop_params.beta = [0.8, 0.7, 0.6];

    variants = {
        'proposed', 'Model 0 (Proposed) - 完全体 (Multi-Subpop + E-Levy + Adaptive + MO-Stop)';
        'no_subpop', 'Model 1 (w/o Multi-Subpop) - 去掉三子种群';
        'no_adaptive', 'Model 2 (w/o Adaptive) - 去掉自适应参数';
        'no_levy', 'Model 3 (w/o E-Levy) - 去掉能量感知Levy跃迁';
        'no_stop', 'Model 4 (w/o MO-SmartStop) - 去掉多目标智能停止'
    };

    reference_point = [1.0, 100000];

    results = struct();
    results.map_name = map_name;
    results.N_User = N_User;
    results.N_UAV = N_UAV;

    for v_idx = 1:size(variants, 1)
        variant_name = variants{v_idx, 1};
        variant_desc = variants{v_idx, 2};
        fprintf('\n--- 运行变体: %s ---\n', variant_desc);

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
        temp_pareto_fronts = cell(n_runs, 1);

        if isempty(gcp('nocreate'))
            fprintf('启动并行池 (%d workers)...\n', n_workers);
            parpool('local', n_workers);
        end

        parfor run = 1:n_runs
            fprintf('  Worker 正在处理 Run %d/%d...\n', run, n_runs);

            stream = RandStream('mt19937ar', 'Seed', run * 100 + v_idx * 1000);
            RandStream.setGlobalStream(stream);

            best_fit = 0;
            bestUAV = zeros(N_UAV, 2);
            cg_curve = zeros(1, 300);
            pareto_archive = struct('Coverage', {}, 'Energy', {}, 'UAV_pos', {});

            switch variant_name
                case 'proposed'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        model0_proposed(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'no_subpop'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        model2_no_subpop(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'no_adaptive'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        model6_no_adaptive(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'no_levy'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        model1_no_levy(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
                case 'no_stop'
                    [best_fit, bestUAV, cg_curve, energy_consumption, pareto_archive] = ...
                        model4_no_stop(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, Ub, Lb, params, priorities);
            end

            best_fits(run) = best_fit;

            center_point = repmat([500, 500], N_UAV, 1);
            fly_dist = sqrt(sum((bestUAV - center_point).^2, 2));
            energies(run) = sum(params.k_move * fly_dist);

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
                temp_pareto_fronts{run} = pareto_front;

                try
                    % HV量纲归一化：能耗目标归一化到[0,1]区间
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
                temp_pareto_fronts{run} = [];
            end

            fprintf('  -> Run %d 完成: Fitness=%.2f, HV=%.4f\n', run, best_fit, hv_values(run));
        end

        pareto_fronts{v_idx} = temp_pareto_fronts;

        results.(variant_name) = struct();
        results.(variant_name).description = variant_desc;
        results.(variant_name).best_fits = best_fits;
        results.(variant_name).energies = energies;
        results.(variant_name).cov_high = cov_high;
        results.(variant_name).cov_total = cov_total;
        results.(variant_name).iter_counts = iter_counts;
        results.(variant_name).convergence_curves = convergence_curves;
        results.(variant_name).mean_fitness = mean(best_fits);
        results.(variant_name).std_fitness = std(best_fits);
        results.(variant_name).mean_energy = mean(energies);
        results.(variant_name).mean_cov_high = mean(cov_high);
        results.(variant_name).mean_cov_total = mean(cov_total);
        results.(variant_name).hv_values = hv_values;
        results.(variant_name).mean_hv = mean(hv_values);
        results.(variant_name).std_hv = std(hv_values);
        results.(variant_name).igd_values = igd_values;
        results.(variant_name).mean_igd = mean(igd_values);
        results.(variant_name).spread_values = spread_values;
        results.(variant_name).mean_spread = mean(spread_values);
        results.(variant_name).pareto_fronts = temp_pareto_fronts;
        results.(variant_name).mean_pareto_size = mean(pareto_sizes);

        fprintf('  >> %s 平均结果:\n', variant_desc);
        fprintf('     平均适应度: %.2f +/- %.2f\n', mean(best_fits), std(best_fits));
        fprintf('     平均能耗: %.2f J\n', mean(energies));
        fprintf('     平均高优先级覆盖率: %.2f%%\n', mean(cov_high));
        fprintf('     平均全局覆盖率: %.2f%%\n', mean(cov_total));
        fprintf('     平均HV: %.4f +/- %.4f\n', mean(hv_values), std(hv_values));
        fprintf('     平均Spread: %.4f\n', mean(spread_values));
    end

    fprintf('\n========== 开始执行消融实验 Min-Max 归一化 ==========\n');
    all_raw_points = [];
    var_keys = fieldnames(results);
    valid_vars = {};
    for i = 1:length(var_keys)
        v = var_keys{i};
        if ~strcmp(v, 'map_name') && ~strcmp(v, 'N_User') && ~strcmp(v, 'N_UAV')
            valid_vars{end+1} = v;
            for run = 1:length(results.(v).pareto_fronts)
                pf = results.(v).pareto_fronts{run};
                if ~isempty(pf)
                    all_raw_points = [all_raw_points; pf];
                end
            end
        end
    end

    c_max = max(all_raw_points(:, 1)); c_min = min(all_raw_points(:, 1));
    e_max = max(all_raw_points(:, 2)); e_min = min(all_raw_points(:, 2));
    if c_max == c_min, c_max = c_min + 1e-6; end
    if e_max == e_min, e_max = e_min + 1e-6; end

    all_points_norm = [];
    for i = 1:length(valid_vars)
        v = valid_vars{i};
        for run = 1:length(results.(v).pareto_fronts)
            pf = results.(v).pareto_fronts{run};
            if ~isempty(pf)
                norm_c = (c_max - pf(:, 1)) / (c_max - c_min);
                norm_e = (pf(:, 2) - e_min) / (e_max - e_min);
                pf_norm = [norm_c, norm_e];
                all_points_norm = [all_points_norm; pf_norm];
                results.(v).pareto_fronts_norm{run} = pf_norm;
            else
                results.(v).pareto_fronts_norm{run} = [];
            end
        end
    end

    all_points_norm = unique(all_points_norm, 'rows');
    true_front_norm = extractNonDominated(all_points_norm);
    ref_point_norm = [1.05, 1.05];

    fprintf('\n========== 消融实验最终汇总表格 (Min-Max 归一化) ==========\n');
    fprintf('%-12s | %-10s | %-10s | %-13s | %-13s | %-8s\n', ...
        'Variant', 'Fitness', 'Energy(J)', 'HV(Norm)↑', 'IGD(Norm)↓', 'Spread↑');
    fprintf('%s\n', repmat('-', 1, 95));

    for i = 1:length(valid_vars)
        v = valid_vars{i};
        r = results.(v);
        hvs = zeros(n_runs, 1); igds = zeros(n_runs, 1); spreads = zeros(n_runs, 1);
        for run = 1:length(r.pareto_fronts_norm)
            pf_norm = r.pareto_fronts_norm{run};
            if ~isempty(pf_norm)
                try
                    metrics = calculate_all_metrics(pf_norm, true_front_norm, ref_point_norm);
                    hvs(run) = metrics.hv; igds(run) = metrics.igd; spreads(run) = metrics.spread;
                catch
                    hvs(run) = NaN; igds(run) = NaN; spreads(run) = NaN;
                end
            end
        end
        results.(v).mean_hv_norm = nanmean(hvs); results.(v).mean_igd_norm = nanmean(igds); results.(v).mean_spread_norm = nanmean(spreads);

        fprintf('%-12s | %-10.2f | %-10.2f | %-5.4f±%-5.4f | %-5.4f±%-5.4f | %-8.4f\n', ...
            v, r.mean_fitness, r.mean_energy, nanmean(hvs), nanstd(hvs), nanmean(igds), nanstd(igds), nanmean(spreads));
    end
    fprintf('%s\n', repmat('-', 1, 95));

    results_file = fullfile(project_dir, 'experiments', ['ablation_results_para_', map_name, '_', datestr(now, 'yyyymmdd_HHMMSS'), '.mat']);
    save(results_file, 'results');
    fprintf('\n完美消融结果已保存至: %s\n', results_file);
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

function non_dominated = extractNonDominated(points)
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

    non_dominated = points(~is_dominated, :);
end