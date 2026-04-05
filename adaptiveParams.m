% adaptiveParams.m - 参数自适应机制
% 根据算法运行状态自动调整参数，减少手动设置
%
% 功能：
% 1. 基于算法性能自动调整参数
% 2. 减少需要手动设置的超参数数量
% 3. 提升算法鲁棒性

function adjusted_params = adaptiveParams(current_params, algorithm_state, iter, FES_max)
    % 根据算法状态自适应调整参数
    %
    % 输入：
    %   current_params: 当前参数结构体
    %   algorithm_state: 算法运行状态（包含适应度、多样性等信息）
    %   iter: 当前迭代次数
    %   FES_max: 最大迭代次数
    %
    % 输出：
    %   adjusted_params: 调整后的参数
    
    adjusted_params = current_params;
    progress = iter / FES_max;

    % ========================================================
    % 修复：注释掉G_weights的自适应调整，避免"伪进化"
    % 绝对不能在运行中途改变目标函数的权重，否则会导致：
    % 1. 相同物理位置因权重变化产生不同分数
    % 2. 误判stagnation_generations，导致精准打击机制失效
    % 3. 适应度虚高但无人机实际未移动
    %
    % 1. 自适应调整子种群权重（基于各子种群性能）
    % if isfield(algorithm_state, 'subpop_fits') && length(algorithm_state.subpop_fits) == 3
    %     subpop_fits = algorithm_state.subpop_fits;
    %     total_fit = sum(subpop_fits);
    %
    %     if total_fit > 0
    %         performance_weights = subpop_fits(:)' / total_fit;
    %
    %         if isfield(current_params, 'G_weights')
    %             alpha = 0.3;
    %             current_weights = current_params.G_weights(:)';
    %             adjusted_params.G_weights = (1 - alpha) * current_weights + ...
    %                 alpha * performance_weights;
    %             adjusted_params.G_weights = adjusted_params.G_weights / sum(adjusted_params.G_weights);
    %             adjusted_params.G_weights = adjusted_params.G_weights(:)';
    %         end
    %     end
    % end
    % ========================================================
    
    % 2. 自适应调整候选解数量K（基于收敛状态）
    if isfield(algorithm_state, 'convergence_rate')
        convergence_rate = algorithm_state.convergence_rate;
        
        % 如果收敛缓慢，增加候选解数量
        if convergence_rate < 0.01 && progress < 0.7
            if isfield(current_params, 'K')
                adjusted_params.K = round(min(current_params.K * 1.1, 60));  % 上限60，确保整数
            end
        % 如果收敛过快，减少候选解数量（加快速度）
        elseif convergence_rate > 0.05 && progress > 0.5
            if isfield(current_params, 'K')
                adjusted_params.K = round(max(current_params.K * 0.95, 20));  % 下限20，确保整数
            end
        end
    end
    
    % 3. 自适应调整U型俯冲概率q（基于探索-开发平衡）
    if isfield(algorithm_state, 'diversity')
        diversity = algorithm_state.diversity;
        
        % 多样性低时，增加探索（降低q，增加V型俯冲）
        if diversity < 0.1 && progress < 0.6
            if isfield(current_params, 'subpop_params') && ...
                    isfield(current_params.subpop_params, 'q')
                adjusted_params.subpop_params.q = ...
                    max(current_params.subpop_params.q * 0.95, 0.2);
            end
        % 多样性高时，增加开发（提高q，增加U型俯冲）
        elseif diversity > 0.3 && progress > 0.4
            if isfield(current_params, 'subpop_params') && ...
                    isfield(current_params.subpop_params, 'q')
                adjusted_params.subpop_params.q = ...
                    min(current_params.subpop_params.q * 1.05, 0.8);
            end
        end
    end
    
    % 4. 自适应调整捕获能力阈值c（基于停滞状态）
    if isfield(algorithm_state, 'stagnation_count')
        stagnation = algorithm_state.stagnation_count;
        
        % 停滞时，降低阈值，增加Levy飞行（探索）
        if stagnation > 15
            if isfield(current_params, 'subpop_params') && ...
                    isfield(current_params.subpop_params, 'c')
                adjusted_params.subpop_params.c = ...
                    max(current_params.subpop_params.c * 0.9, 0.1);
            end
        % 快速改进时，提高阈值，增加转向（开发）
        elseif stagnation < 5 && progress > 0.3
            if isfield(current_params, 'subpop_params') && ...
                    isfield(current_params.subpop_params, 'c')
                adjusted_params.subpop_params.c = ...
                    min(current_params.subpop_params.c * 1.1, 0.5);
            end
        end
    end
    
    % 5. 自适应调整初始标准差sigma0（基于搜索范围）
    if isfield(algorithm_state, 'search_range')
        search_range = algorithm_state.search_range;
        
        % 搜索范围小，增加sigma0（扩大搜索）
        if search_range < 0.1 * norm([1000, 1000]) && progress < 0.5
            if isfield(current_params, 'subpop_params') && ...
                    isfield(current_params.subpop_params, 'sigma0')
                adjusted_params.subpop_params.sigma0 = ...
                    min(current_params.subpop_params.sigma0 * 1.1, 50);
            end
        % 搜索范围大，减少sigma0（聚焦搜索）
        elseif search_range > 0.3 * norm([1000, 1000]) && progress > 0.6
            if isfield(current_params, 'subpop_params') && ...
                    isfield(current_params.subpop_params, 'sigma0')
                adjusted_params.subpop_params.sigma0 = ...
                    max(current_params.subpop_params.sigma0 * 0.95, 10);
            end
        end
    end
    
    % 6. 基于迭代进度的自适应调整
    % 前期：强调探索
    if progress < 0.3
        if isfield(current_params, 'subpop_params')
            % 降低q（更多V型俯冲）
            if isfield(current_params.subpop_params, 'q')
                adjusted_params.subpop_params.q = ...
                    current_params.subpop_params.q * 0.9;
            end
            % 增加sigma0（更大搜索范围）
            if isfield(current_params.subpop_params, 'sigma0')
                adjusted_params.subpop_params.sigma0 = ...
                    current_params.subpop_params.sigma0 * 1.1;
            end
        end
    % 后期：强调开发
    elseif progress > 0.7
        if isfield(current_params, 'subpop_params')
            % 提高q（更多U型俯冲）
            if isfield(current_params.subpop_params, 'q')
                adjusted_params.subpop_params.q = ...
                    min(current_params.subpop_params.q * 1.1, 0.8);
            end
            % 减少sigma0（聚焦搜索）
            if isfield(current_params.subpop_params, 'sigma0')
                adjusted_params.subpop_params.sigma0 = ...
                    max(current_params.subpop_params.sigma0 * 0.9, 10);
            end
        end
    end
end

% 计算算法状态指标
function state = computeAlgorithmState(mem_matrix, prev_fits, curr_fits, ...
    best_fit, iter, FES_max, User, priorities, params)
    % 计算算法运行状态指标
    
    state = struct();
    
    % 1. 收敛速度（适应度改进率）
    if length(prev_fits) > 0 && length(curr_fits) > 0
        prev_best = max(prev_fits);
        curr_best = max(curr_fits);
        if prev_best > 0
            state.convergence_rate = abs(curr_best - prev_best) / prev_best;
        else
            state.convergence_rate = 0;
        end
    else
        state.convergence_rate = 0;
    end
    
    % 2. 多样性（基于适应度方差）
    all_fits = [];
    for g = 1:length(mem_matrix)
        for i = 1:size(mem_matrix{g}, 1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1 && size(candidate, 2) == size(mem_matrix{g}, 3) * 2
                candidate = reshape(candidate, size(mem_matrix{g}, 3), 2);
            end
            % 简化：使用位置方差作为多样性指标
            if size(candidate, 1) > 1
                pos_variance = var(candidate(:));
                all_fits = [all_fits, pos_variance];
            end
        end
    end
    if length(all_fits) > 1
        state.diversity = std(all_fits) / (mean(all_fits) + 1e-6);
    else
        state.diversity = 0.5;
    end
    
    % 3. 停滞计数（连续无改进代数）
    state.stagnation_count = 0;  % 需要外部传入
    
    % 4. 适应度改进
    if length(prev_fits) > 0 && length(curr_fits) > 0
        state.fitness_improvement = max(curr_fits) - max(prev_fits);
    else
        state.fitness_improvement = 0;
    end
    
    % 5. 搜索范围（UAV位置分布范围）
    all_positions = [];
    for g = 1:length(mem_matrix)
        for i = 1:size(mem_matrix{g}, 1)
            candidate = squeeze(mem_matrix{g}(i,:,:));
            if size(candidate, 1) == 1
                candidate = reshape(candidate, size(mem_matrix{g}, 3), 2);
            end
            all_positions = [all_positions; candidate];
        end
    end
    if size(all_positions, 1) > 1
        pos_range = max(all_positions) - min(all_positions);
        state.search_range = norm(pos_range);
    else
        state.search_range = 0;
    end
    
    % 6. 子种群适应度
    state.subpop_fits = curr_fits;
end

