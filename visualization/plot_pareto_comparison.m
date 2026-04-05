function plot_pareto_comparison(data_file, output_dir)
if nargin < 1
    data_file = fullfile('..', 'experiments', 'comparison_results_Map1_Medium_20260403_173104.mat');
end
if nargin < 2
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4]);
set(gcf, 'Color', 'w');

%% Figure A: Comparative Analysis - All algorithms
ax1 = subplot(1, 2, 1);

algo_names = {'PSO', 'GA', 'GOA', 'GWO', 'cSA', 'Proposed'};
algo_colors = {[0.3, 0.3, 0.3], [0.5, 0.5, 0.5], [0.6, 0.6, 0.6], [0.7, 0.7, 0.7], [0.8, 0.8, 0.8], [1, 0, 0]};
algo_markers = {'o', 's', '^', 'd', 'v', '*'};

if exist('results', 'var')
    algos = fieldnames(results);
    for i = 1:min(length(algos), 6)
        algo = algos{i};
        pf_data = results.(algo).pareto_fronts;
        if ~isempty(pf_data) && length(pf_data) >= 1 && ~isempty(pf_data{1})
            all_pf = [];
            for run = 1:min(length(pf_data), 10)
                if ~isempty(pf_data{run})
                    all_pf = [all_pf; pf_data{run}];
                end
            end
            if ~isempty(all_pf)
                hold on;
                if i == 6
                    scatter(all_pf(:,1)*100, all_pf(:,2)/1000, 80, algo_colors{min(i,6)}, algo_markers{min(i,6)}, 'filled', 'LineWidth', 1.5);
                else
                    scatter(all_pf(:,1)*100, all_pf(:,2)/1000, 30, algo_colors{min(i,6)}, algo_markers{min(i,6)}, 'filled');
                end
            end
        end
    end
elseif exist('pareto_fronts', 'var')
    for i = 1:min(length(pareto_fronts), 6)
        pf = pareto_fronts{i};
        if ~isempty(pf)
            hold on;
            if i == 6
                scatter(pf(:,1)*100, pf(:,2)/1000, 80, algo_colors{i}, algo_markers{i}, 'filled', 'LineWidth', 1.5);
            else
                scatter(pf(:,1)*100, pf(:,2)/1000, 30, algo_colors{i}, algo_markers{i}, 'filled');
            end
        end
    end
end

xlabel('Coverage (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy Consumption (kJ)', 'FontSize', 12, 'FontWeight', 'bold');
title('(a) Pareto Front Comparison', 'FontSize', 14, 'FontWeight', 'bold');
legend(algo_names{1:6}, 'Location', 'best', 'FontSize', 9);
grid on;
set(gca, 'FontSize', 10);

%% Figure B: Ablation Study - Proposed vs no_subpop
ax2 = subplot(1, 2, 2);

ablation_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_20260405_184023.mat');
if exist(ablation_file, 'file')
    load(ablation_file);
else
    ablation_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_20260405_001223.mat');
    if exist(ablation_file, 'file')
        load(ablation_file);
    end
end

if exist('results', 'var')
    variants = {'proposed', 'no_subpop'};
    variant_labels = {'Proposed (Multi-Subpop)', 'w/o Multi-Subpop'};
    variant_colors = {[1, 0, 0], [0.5, 0.5, 0.5]};
    variant_markers = {'*', 'o'};

    for v = 1:length(variants)
        variant = variants{v};
        if isfield(results, variant)
            pf_data = results.(variant).pareto_fronts;
            if ~isempty(pf_data) && length(pf_data) >= 1
                all_pf = [];
                for run = 1:min(length(pf_data), 10)
                    if ~isempty(pf_data{run})
                        all_pf = [all_pf; pf_data{run}];
                    end
                end
                if ~isempty(all_pf)
                    hold on;
                    if v == 1
                        scatter(all_pf(:,1)*100, all_pf(:,2)/1000, 100, variant_colors{v}, variant_markers{v}, 'filled', 'LineWidth', 1.5);
                    else
                        scatter(all_pf(:,1)*100, all_pf(:,2)/1000, 30, variant_colors{v}, variant_markers{v}, 'filled');
                    end
                end
            end
        end
    end

    legend(variant_labels, 'Location', 'best', 'FontSize', 10);
end

xlabel('Coverage (%)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Energy Consumption (kJ)', 'FontSize', 12, 'FontWeight', 'bold');
title('(b) Ablation Study: Multi-Subpop Effect', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 10);

annotation('textarrow', [0.35, 0.25], [0.15, 0.1], 'String', 'Subpop divides & conquers', 'FontSize', 10, 'Color', [0.3, 0.3, 0.3]);

saveas(fig, fullfile(output_dir, 'pareto_comparison.fig'));
saveas(fig, fullfile(output_dir, 'pareto_comparison.png'));
fprintf('Pareto comparison saved to %s\n', output_dir);
close(fig);
end