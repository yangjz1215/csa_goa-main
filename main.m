% main.m - 单目标优化主程序
% 优化目标：最大化适应度（优先级覆盖之和）

clear; clc; close all;
rng('shuffle');  % 根据当前时间生成真随机数

%% 1. 场景参数配置
Lb = [0, 0];
Ub = [1000, 1000];

N_User = 500;
User = Lb + (Ub - Lb) .* rand(N_User, 2);

ratios = [0.25, 0.25, 0.25, 0.25];
num_levels = length(ratios);
counts = floor(ratios * N_User);
remainder = N_User - sum(counts);
extra_indices = randperm(num_levels, remainder);
for i = 1:remainder
    counts(extra_indices(i)) = counts(extra_indices(i)) + 1;
end
task_levels = [];
for i = 1:num_levels
    task_levels = [task_levels; i * ones(counts(i), 1)];
end
priorities = task_levels(randperm(N_User));
U_high_idx = find(priorities >= 3);

%% 1.5 创新点：长尾用户机会成本评估 (孤岛衰减机制)
% 计算每个用户的有效优先级 (Effective Priority)
center_point = [500, 500];
cover_radius = 150;  % 覆盖半径（与params保持一致）
effective_priorities = zeros(N_User, 1);

for i = 1:N_User
    dists_to_others = sqrt(sum((User - User(i,:)).^2, 2));
    neighbor_count = sum(dists_to_others <= cover_radius) - 1;

    dist_to_center = norm(User(i,:) - center_point);

    if neighbor_count == 0
        discount_factor = exp(-dist_to_center / 500);
    elseif neighbor_count < 3
        discount_factor = 0.8 + 0.2 * exp(-dist_to_center / 800);
    else
        discount_factor = 1.0;
    end

    effective_priorities(i) = priorities(i) * discount_factor;
end

priorities = effective_priorities;
U_high_idx = find(priorities >= 3);

N_RRH = 10;
N_eRRH = 4;
RRH = GenerateRRH(N_RRH, Ub, Lb);
RRH_type = zeros(N_RRH, 1);
RRH_type(1:N_eRRH) = 1;

N_UAV = 15;
N_eUAV = 5;
UAV_type = zeros(N_UAV, 1);
UAV_type(1:N_eUAV) = 1;

params = struct();
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

D_max = 500;
D_min = 100;
params.D = ((D_max-D_min)*rand(N_User,1)+D_min)*8192;

C_max = 1.0;
C_min = 0.5;
params.C = ((C_max-C_min)*rand(N_User,1)+C_min)*10^9;

DT_max = 1.2;
DT_min = 0.8;
params.DT = ((DT_max-DT_min)*rand(N_User,1)+DT_min);

%% 2. 算法参数配置
% ========== 超参数优化方法选择 ==========
% 可选方法：
%   'fixed'             - 固定参数（手动设置所有参数）
%   'adaptive'          - 自适应调整（推荐日常使用，自动调整参数）
%   'bayesian'          - 贝叶斯优化（离线优化，找到最优参数，耗时较长）
%   'bayesian_adaptive' - 贝叶斯优化初始化 + 自适应调整（推荐！兼顾初始质量和实时调整）
optimization_method = 'bayesian_adaptive';

params.FES_max = 300;
params.K = 40;

switch optimization_method
    case 'fixed'
        params.G_weights = [0.35, 0.45, 0.2];
        params.subpop_params = struct(...
            'q', [0.5, 0.4, 0.3], ...
            'c', [0.3, 0.25, 0.2], ...
            'beta', [0.2, 0.3, 0.1], ...
            'sigma0', [25, 30, 20]);
        params.enable_adaptive = false;
        fprintf('使用固定参数方法\n');
        
    case 'adaptive'
        params.G_weights = [0.35, 0.45, 0.2];
        params.subpop_params = struct(...
            'q', [0.5, 0.4, 0.3], ...
            'c', [0.3, 0.25, 0.2], ...
            'beta', [0.2, 0.3, 0.1], ...
            'sigma0', [25, 30, 20]);
        params.enable_adaptive = true;
        fprintf('使用自适应参数调整方法\n');
        
    case 'bayesian'
        fprintf('========== 使用贝叶斯优化 ==========\n');
        fprintf('这将运行多次算法评估以找到最优参数，可能需要较长时间...\n');
        
        scenario_config = struct();
        scenario_config.N_User = N_User;
        scenario_config.User = User;
        scenario_config.N_RRH = N_RRH;
        scenario_config.RRH = RRH;
        scenario_config.RRH_type = RRH_type;
        scenario_config.N_UAV = N_UAV;
        scenario_config.UAV_type = UAV_type;
        scenario_config.Ub = Ub;
        scenario_config.Lb = Lb;
        scenario_config.priorities = priorities;
        
        options = struct();
        options.max_iterations = 15;
        options.initial_samples = 5;
        
        best_params = optimizeHyperparameters('bayesian', scenario_config, options);
        
        params.FES_max = best_params.FES_max;
        params.K = best_params.K;
        params.G_weights = best_params.G_weights;
        params.subpop_params = best_params.subpop_params;
        params.enable_adaptive = false;
        
    case 'bayesian_adaptive'
        fprintf('========== 使用贝叶斯优化初始化 + 自适应调整 ==========\n');
        
        bayesian_params_file = 'bayesian_optimized_params.mat';
        
        if exist(bayesian_params_file, 'file')
            fprintf('加载已保存的贝叶斯优化参数...\n');
            loaded = load(bayesian_params_file);
            best_params = loaded.best_params;
        else
            fprintf('首次运行，执行贝叶斯优化...\n');
            fprintf('这可能需要较长时间，优化完成后参数将保存以供后续使用\n');
            
            scenario_config = struct();
            scenario_config.N_User = N_User;
            scenario_config.User = User;
            scenario_config.N_RRH = N_RRH;
            scenario_config.RRH = RRH;
            scenario_config.RRH_type = RRH_type;
            scenario_config.N_UAV = N_UAV;
            scenario_config.UAV_type = UAV_type;
            scenario_config.Ub = Ub;
            scenario_config.Lb = Lb;
            scenario_config.priorities = priorities;
            
            options = struct();
            options.max_iterations = 15;
            options.initial_samples = 5;
            
            best_params = optimizeHyperparameters('bayesian', scenario_config, options);
            
            save(bayesian_params_file, 'best_params');
            fprintf('贝叶斯优化参数已保存到 %s\n', bayesian_params_file);
        end
        
        params.FES_max = best_params.FES_max;
        params.K = best_params.K;
        params.G_weights = best_params.G_weights;
        params.subpop_params = best_params.subpop_params;
        params.enable_adaptive = true;
        
        fprintf('贝叶斯优化初始化完成，启用自适应调整\n');
        
    otherwise
        error('未知的优化方法: %s。请选择: fixed, adaptive, bayesian, bayesian_adaptive', optimization_method);
end

%% 3. 运行算法
fprintf('\n开始运行cSA-GOA算法（方法: %s）...\n', optimization_method);
[best_fit, bestUAV, cg_curve, energy_consumption, E_remaining_history, final_E_remaining, curr_curve, actual_iter, ~, pareto_archive] = ...
    cSA_GOA_main(N_User, User, N_RRH, RRH, RRH_type, N_UAV, UAV_type, ...
    Ub, Lb, params, priorities);

%% 4. 结果可视化
figure('Name','收敛曲线','Position',[100,100,800,500]);
plot(1:actual_iter, curr_curve(1:actual_iter), 'b-','LineWidth',1.5);
xlabel('迭代次数'); ylabel('当前适应度');
title('算法收敛曲线'); 
grid on;
xticks(0:10:actual_iter);
xlim([1, actual_iter]);

figure('Name','能耗变化','Position',[200,200,800,500]);
plot(1:actual_iter, energy_consumption(1:actual_iter), 'r-','LineWidth',1.5);
xlabel('迭代次数'); ylabel('总能耗（J）');
title('无人机总能耗变化'); grid on;

% 计算能耗用于Pareto图
center_point = repmat([500, 500], N_UAV, 1);
fly_distances = sqrt(sum((bestUAV - center_point).^2, 2));
total_energy_consumed = params.k_move * sum(fly_distances);
final_cov_total = calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius);

% Pareto前沿分析
if ~isempty(pareto_archive)
    arch_cov = [pareto_archive.Coverage];
    arch_energy = [pareto_archive.Energy];

    figure('Name','Pareto 前沿分析','Position',[400,200,600,500]);
    scatter(arch_cov, arch_energy, 60, 'filled', 'MarkerFaceColor',[0.8500 0.3250 0.0980]);
    hold on;

    [sorted_cov, sort_idx] = sort(arch_cov);
    sorted_energy = arch_energy(sort_idx);
    plot(sorted_cov, sorted_energy, 'b-', 'LineWidth', 1.5);

    grid on;
    title('无人机部署的多目标 Pareto 前沿 (被动归档)');
    xlabel('总覆盖用户数 (越大越好)');
    ylabel('系统总飞行能耗 / J (越小越好)');
    scatter(final_cov_total * N_User, total_energy_consumed, 150, 'p', 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k');

    % ==========================================================
    % 论文亮点：后验多目标决策分析 (A Posteriori Decision Making)
    % 从被动收集的 Pareto 前沿中提取代表性决策方案
    % ==========================================================

    [max_cov, idx_max_cov] = max(arch_cov);
    [min_eng, idx_min_eng] = min(arch_energy);

    solution_A = pareto_archive(idx_max_cov);
    solution_B = pareto_archive(idx_min_eng);

    norm_cov = (arch_cov - min(arch_cov)) / (max(arch_cov) - min(arch_cov) + 1e-6);
    norm_eng = (arch_energy - min(arch_energy)) / (max(arch_energy) - min(arch_energy) + 1e-6);
    distances_to_ideal = sqrt((1 - norm_cov).^2 + (0 - norm_eng).^2);
    [~, idx_knee] = min(distances_to_ideal);
    solution_Knee = pareto_archive(idx_knee);

    figure(gcf);
    hold on;
    scatter(solution_Knee.Coverage, solution_Knee.Energy, 200, 'pentagram', 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    scatter(solution_A.Coverage, solution_A.Energy, 150, 'd', 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k');
    scatter(solution_B.Coverage, solution_B.Energy, 150, 'd', 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'k');

    legend('Pareto非支配解', 'Pareto前沿边界', '单目标加权最终解', ...
           'Knee Point (最高性价比解)', '最大覆盖解 (激进)', '最低能耗解 (保守)', 'Location', 'best');

    fprintf('\n========== 多目标决策方案推荐 (A Posteriori) ==========\n');
    fprintf('【方案0：原单目标加权解 (黄色星星)】\n');
    fprintf('  - 总覆盖人数: %d, 总能耗: %.2f J\n', final_cov_total * N_User, total_energy_consumed);
    fprintf('  - 特点: 照顾高优先级与分散度，无视物理前沿性价比\n\n');

    fprintf('【方案1：Knee Point 最高性价比解 (绿色星星)】\n');
    fprintf('  - 总覆盖人数: %d, 总能耗: %.2f J\n', solution_Knee.Coverage, solution_Knee.Energy);
    fprintf('  - 特点: Pareto 前沿上距理想目标最近的点，权衡最佳\n\n');

    fprintf('【方案2：最大覆盖激进解 (品红菱形)】\n');
    fprintf('  - 总覆盖人数: %d, 总能耗: %.2f J\n', solution_A.Coverage, solution_A.Energy);
    fprintf('  - 特点: 牺牲能耗，榨干最后一点覆盖率\n\n');

    fprintf('【方案3：最低能耗保守解 (青色菱形)】\n');
    fprintf('  - 总覆盖人数: %d, 总能耗: %.2f J\n', solution_B.Coverage, solution_B.Energy);
    fprintf('  - 特点: 能耗最低，覆盖人数也相当可观\n\n');
    fprintf('=========================================================\n');
end

figure('Name','最终部署','Position',[300,300,900,700]);
plot(User(priorities<3,1), User(priorities<3,2), 'bo','MarkerSize',6,'DisplayName','低优先级用户');
hold on;
plot(User(priorities>=3,1), User(priorities>=3,2), 'ro','MarkerSize',8,'DisplayName','高优先级用户');
plot(RRH(:,1), RRH(:,2), 'g^','MarkerSize',10,'DisplayName','RRH');
plot(bestUAV(UAV_type==0,1), bestUAV(UAV_type==0,2), 'k*','MarkerSize',12,'DisplayName','普通UAV');
plot(bestUAV(UAV_type==1,1), bestUAV(UAV_type==1,2), 'm*','MarkerSize',12,'DisplayName','增强UAV');
xlabel('X坐标（米）'); ylabel('Y坐标（米）');
legend('Location','best'); axis equal; grid on;

fprintf('\n========== 最终优化结果汇总 ==========\n');
fprintf('优化方法: %s\n', optimization_method);
fprintf('1. 历史最优适应度：%.4f\n', best_fit);
fprintf('2. 最终当前适应度：%.4f\n', curr_curve(actual_iter));
fprintf('3. 最终总能耗：%.2f J（剩余：%.2f J）\n', energy_consumption(actual_iter), sum(final_E_remaining));
fprintf('4. 高优先级覆盖率：%.2f%%\n', calcCoverageWithRRH(bestUAV, User(U_high_idx,:), params.cover_radius, RRH, params.RRH_radius)*100);
fprintf('5. 全局覆盖率：%.2f%%\n', calcCoverageWithRRH(bestUAV, User, params.cover_radius, RRH, params.RRH_radius)*100);
fprintf('6. Pareto归档非支配解数量：%d\n', length(pareto_archive));
if ~isempty(pareto_archive)
    pareto_cov = [pareto_archive.Coverage];
    pareto_energy = [pareto_archive.Energy];
    fprintf('   最优覆盖方案：%d人（能耗%.2f J）\n', max(pareto_cov), pareto_energy(pareto_cov==max(pareto_cov)));
    fprintf('   最节能方案：%.2f J（覆盖%d人）\n', min(pareto_energy), pareto_cov(pareto_energy==min(pareto_energy)));
end
fprintf('=====================================\n');

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
