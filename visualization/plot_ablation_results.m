function plot_ablation_results(data_file, output_dir)
if nargin < 1
    data_file = fullfile('..', 'experiments', 'ablation_results_para_Map1_Medium_20260405_184023.mat');
end
if nargin < 2
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.45]);
set(gcf, 'Color', 'w');

variants = {'proposed', 'no_subpop', 'no_adaptive', 'no_levy', 'no_stop'};
variant_labels = {'Proposed', 'w/o Multi-Subpop', 'w/o Adaptive', 'w/o E-Levy', 'w/o MO-SmartStop'};

hv_means = zeros(length(variants), 1);
hv_stds = zeros(length(variants), 1);
iter_means = zeros(length(variants), 1);
iter_stds = zeros(length(variants), 1);

for v = 1:length(variants)
    variant = variants{v};
    if isfield(results, variant)
        hv_means(v) = results.(variant).mean_hv;
        hv_stds(v) = results.(variant).std_hv;
        iter_means(v) = mean(results.(variant).iter_counts);
        iter_stds(v) = std(results.(variant).iter_counts);
    end
end

x = 1:length(variants);
width = 0.35;

hv_normalized = hv_means / max(hv_means) * 100;
iter_normalized = (1 - iter_means / max(iter_means)) * 100;

bar1 = bar(x - width/2, hv_normalized, width, 'FaceColor', [0.85, 0.33, 0.31], 'EdgeColor', [0.6, 0.2, 0.1], 'LineWidth', 1.5);
hold on;
bar2 = bar(x + width/2, iter_normalized, width, 'FaceColor', [0.33, 0.59, 0.74], 'EdgeColor', [0.2, 0.4, 0.6], 'LineWidth', 1.5);

for i = 1:length(variants)
    text(x(i) - width/2, hv_normalized(i) + 2, sprintf('%.3f', hv_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.6, 0.2, 0.1]);
    text(x(i) + width/2, iter_normalized(i) + 2, sprintf('%.0f', iter_means(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.2, 0.4, 0.6]);
end

set(gca, 'XTick', x);
set(gca, 'XTickLabel', variant_labels, 'FontSize', 9);
xtickangle(45);

ylabel('Normalized Score (%)', 'FontSize', 12, 'FontWeight', 'bold');

title('Ablation Study: Performance vs Efficiency Trade-off', 'FontSize', 14, 'FontWeight', 'bold');

legend([bar1, bar2], {'Hypervolume (higher \uparrow)', 'Efficiency (lower iter \downarrow)'}, ...
    'Location', 'northoutside', 'FontSize', 10, 'Orientation', 'horizontal');

grid on;
set(gca, 'YGrid', 'on');
set(gca, 'YMinorGrid', 'off');
ylim([0, 120]);

annotation('textbox', [0.15, 0.25, 0.2, 0.1], 'String', {'Proposed achieves', 'near-optimal HV with', 'significantly fewer iterations'}, ...
    'Color', [0.85, 0.33, 0.31], 'FontSize', 9, 'FontWeight', 'bold', ...
    'EdgeColor', [0.85, 0.33, 0.31], 'FaceAlpha', 0.1, 'BackgroundColor', [1, 1, 1]);

saveas(fig, fullfile(output_dir, 'ablation_results.fig'));
saveas(fig, fullfile(output_dir, 'ablation_results.png'));
fprintf('Ablation results saved to %s\n', output_dir);
close(fig);

fig2 = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.7, 0.35]);
set(gcf, 'Color', 'w');

subplot(1, 2, 1);
bar(x, hv_means, 'FaceColor', [0.85, 0.33, 0.31], 'EdgeColor', [0.6, 0.2, 0.1], 'LineWidth', 1.5);
hold on;
errorbar(x, hv_means, hv_stds, 'k.', 'LineWidth', 1.5, 'CapSize', 5);
set(gca, 'XTick', x);
set(gca, 'XTickLabel', variant_labels, 'FontSize', 8);
xtickangle(45);
ylabel('Hypervolume', 'FontSize', 11, 'FontWeight', 'bold');
title('HV Comparison (Mean \pm Std)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

subplot(1, 2, 2);
bar(x, iter_means, 'FaceColor', [0.33, 0.59, 0.74], 'EdgeColor', [0.2, 0.4, 0.6], 'LineWidth', 1.5);
hold on;
errorbar(x, iter_means, iter_stds, 'k.', 'LineWidth', 1.5, 'CapSize', 5);
set(gca, 'XTick', x);
set(gca, 'XTickLabel', variant_labels, 'FontSize', 8);
xtickangle(45);
ylabel('Average Iterations', 'FontSize', 11, 'FontWeight', 'bold');
title('Computational Cost (Mean \pm Std)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

saveas(fig2, fullfile(output_dir, 'ablation_detailed.fig'));
saveas(fig2, fullfile(output_dir, 'ablation_detailed.png'));
fprintf('Detailed ablation results saved to %s\n', output_dir);
close(fig2);
end