function plot_ablation_results(data_file, output_dir)
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_20260405_184023.mat');
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

variants = {'proposed', 'no_subpop', 'no_adaptive', 'no_levy', 'no_stop'};
labels = {'Proposed', 'w/o Subpop', 'w/o Adaptive', 'w/o E-Levy', 'w/o MO-Stop'};

hv_means = zeros(1, 5);
iter_means = zeros(1, 5);

for i = 1:5
    v = variants{i};
    if isfield(results, v)
        hv_means(i) = results.(v).mean_hv;
        iter_means(i) = mean(results.(v).iter_counts);
    end
end

fig = figure('Position', [100, 100, 800, 500]);
set(gcf, 'Color', 'w');

yyaxis left
b1 = bar(1:5, hv_means, 0.4, 'FaceColor', [0.0000, 0.4470, 0.7410], 'EdgeColor', 'k', 'LineWidth', 1.2);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold', 'Color', [0.0000, 0.4470, 0.7410]);
ylim([min(hv_means)*0.95, max(hv_means)*1.05]);
set(gca, 'ycolor', [0.0000, 0.4470, 0.7410]);

yyaxis right
hold on;
b2 = bar((1:5)+0.4, iter_means, 0.4, 'FaceColor', [0.8500, 0.3250, 0.0980], 'EdgeColor', 'k', 'LineWidth', 1.2);
ylabel('Actual Iterations (Cost)', 'FontWeight', 'bold', 'Color', [0.8500, 0.3250, 0.0980]);
ylim([0, 350]);
set(gca, 'ycolor', [0.8500, 0.3250, 0.0980]);

xticks((1:5) + 0.2);
xticklabels(labels);

legend([b1, b2], {'Performance (HV)', 'Computational Cost (Iter)'}, 'Location', 'northeast');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on;
box on;

saveas(fig, fullfile(output_dir, 'ablation_detailed.fig'));
saveas(fig, fullfile(output_dir, 'ablation_detailed.png'));
fprintf('Ablation results saved to %s\n', output_dir);
close(fig);
end