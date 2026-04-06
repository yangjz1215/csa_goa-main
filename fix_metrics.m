function fix_metrics(results_file)
    project_dir = fileparts(mfilename('fullpath'));
    if isempty(project_dir)
        project_dir = pwd;
    end

    if nargin < 1
        files = dir(fullfile(project_dir, 'experiments', 'comparison_results_para_Map1_Medium_*.mat'));
        if isempty(files)
            fprintf('ERROR: No comparison results found in experiments/\n');
            return;
        end
        [~, idx] = sort([files.datenum], 'descend');
        results_file = fullfile(project_dir, files(idx(1).folder), files(idx(1).name);
        fprintf('Using latest file: %s\n', results_file);
    else
        if ~contains(results_file, '/') && ~contains(results_file, '\')
            results_file = fullfile(project_dir, 'experiments', results_file);
        end
    end

    load(results_file);

    alg_names = fieldnames(results);
    alg_names = alg_names(~strcmp(alg_names, 'map_name') & ...
                          ~strcmp(alg_names, 'N_User') & ...
                          ~strcmp(alg_names, 'N_UAV'));

    max_energy = 50000 * 15;
    ref_point_norm = [1.0, 1.0];

    fprintf('\n========== 修正量纲后的完美表格 ==========\n');
    fprintf('%-12s | %-10s | %-10s | %-10s | %-8s | %-8s\n', ...
        'Algorithm', 'Fitness', 'Energy(J)', 'HV(Norm)', 'IGD(Norm)', 'Spread');
    fprintf('%s\n', repmat('-', 1, 85));

    for i = 1:length(alg_names)
        alg = alg_names{i};
        r = results.(alg);

        fprintf('%-12s | %-10.2f | %-10.2f | ', ...
            alg, r.mean_fitness, r.mean_energy);

        if isfield(r, 'mean_hv_norm')
            fprintf('%-10.4f | %-10.4f | %-8.4f\n', ...
                r.mean_hv_norm, r.mean_igd_norm, r.mean_spread_norm);
        else
            fprintf('  (需重新计算)\n');
        end
    end

    fprintf('\n========== 正在重新计算归一化指标 ==========\n');

    all_points_norm = [];
    for i = 1:length(alg_names)
        alg = alg_names{i};
        r = results.(alg);

        if isfield(r, 'pareto_fronts')
            for run = 1:length(r.pareto_fronts)
                pf = r.pareto_fronts{run};
                if ~isempty(pf)
                    pf_norm = [pf(:, 1), pf(:, 2) / max_energy];
                    all_points_norm = [all_points_norm; pf_norm];
                end
            end
        end
    end

    if size(all_points_norm, 1) == 0
        fprintf('ERROR: No Pareto front data found!\n');
        fprintf('请确保 run_comparison_para.m 已包含保存 pareto_fronts 的代码，\n');
        fprintf('然后重新运行对比实验生成新的 .mat 文件。\n');
        return;
    end

    true_front_norm = extractNonDominated(all_points_norm);
    fprintf('  归一化后的真实前沿点数: %d\n', size(true_front_norm, 1));

    fprintf('\n========== 修正量纲后的完美表格 ==========\n');
    fprintf('%-12s | %-10s | %-10s | %-12s | %-12s | %-8s\n', ...
        'Algorithm', 'Fitness', 'Energy(J)', 'HV(Norm)↑', 'IGD(Norm)↓', 'Spread↑');
    fprintf('%s\n', repmat('-', 1, 95));

    for i = 1:length(alg_names)
        alg = alg_names{i};
        r = results.(alg);

        hvs = zeros(length(r.pareto_fronts), 1);
        igds = zeros(length(r.pareto_fronts), 1);
        spreads = zeros(length(r.pareto_fronts), 1);

        for run = 1:length(r.pareto_fronts)
            pf = r.pareto_fronts{run};
            if ~isempty(pf)
                pf_norm = [pf(:, 1), pf(:, 2) / max_energy];
                try
                    metrics = calculate_all_metrics(pf_norm, true_front_norm, ref_point_norm);
                    hvs(run) = metrics.hv;
                    igds(run) = metrics.igd;
                    spreads(run) = metrics.spread;
                catch
                    hvs(run) = NaN;
                    igds(run) = NaN;
                    spreads(run) = NaN;
                end
            end
        end

        results.(alg).mean_hv_norm = mean(hvs);
        results.(alg).std_hv_norm = std(hvs);
        results.(alg).mean_igd_norm = mean(igds);
        results.(alg).std_igd_norm = std(igds);
        results.(alg).mean_spread_norm = mean(spreads);
        results.(alg).std_spread_norm = std(spreads);

        fprintf('%-12s | %-10.2f | %-10.2f | %-5.4f±%-5.4f | %-5.4f±%-5.4f | %-8.4f\n', ...
            alg, r.mean_fitness, r.mean_energy, ...
            mean(hvs), std(hvs), mean(igds), std(igds), mean(spreads));
    end

    fprintf('%s\n', repmat('-', 1, 95));
    fprintf('注: HV↑越大越好 | IGD↓越小越好 | Spread↑分布越均匀\n');

    [~, name] = fileparts(results_file);
    save(fullfile(project_dir, 'experiments', [name, '_fixed.mat']), 'results');
    fprintf('\n修正后的结果已保存: %s_fixed.mat\n', name);
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