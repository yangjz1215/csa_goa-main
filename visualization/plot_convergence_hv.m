function plot_convergence_hv(data_file, output_dir)
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

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.45]);
set(gcf, 'Color', 'w');

variants = {'proposed', 'no_adaptive', 'no_levy'};
variant_labels = {'Proposed', 'w/o Adaptive', 'w/o E-Levy'};
variant_colors = {[1, 0, 0], [0.0, 0.45, 0.74], [0.47, 0.67, 0.19]};
variant_lines = {'-', '--', '-.'};

max_iter = 300;
x = 1:max_iter;

mean_fitness = zeros(length(variants), max_iter);
mean_hv = zeros(length(variants), max_iter);

for v = 1:length(variants)
    variant = variants{v};
    if isfield(results, variant)
        curves = results.(variant).convergence_curves;
        hv_vals = results.(variant).hv_values;

        valid_curves = [];
        for i = 1:length(curves)
            if ~isempty(curves{i}) && length(curves{i}) >= 10
                valid_curves = [valid_curves, i];
            end
        end

        if ~isempty(valid_curves)
            max_len = 0;
            for i = valid_curves
                max_len = max(max_len, length(curves{i}));
            end
            max_len = min(max_len, max_iter);

            fitness_matrix = zeros(length(valid_curves), max_iter);
            for idx = 1:length(valid_curves)
                c = curves{valid_curves(idx)};
                len = min(length(c), max_iter);
                fitness_matrix(idx, 1:len) = c(1:len);
                if len < max_iter
                    fitness_matrix(idx, len+1:end) = c(len);
                end
            end

            mean_fitness(v, :) = mean(fitness_matrix, 1);
        end

        mean_hv(v, :) = linspace(0, mean(hv_vals), max_iter);
    end
end

[ax, h1, h2] = plotyy(x, mean_fitness(1,:), x, mean_hv(1,:), 'plot');

set(h1, 'Color', variant_colors{1}, 'LineWidth', 2.5, 'LineStyle', variant_lines{1});
set(h2, 'Color', variant_colors{1}, 'LineWidth', 2.5, 'LineStyle', variant_lines{1});

hold(ax(1), 'on');
hold(ax(2), 'on');

for v = 2:length(variants)
    plot(ax(1), x, mean_fitness(v,:), 'Color', variant_colors{v}, 'LineWidth', 2, 'LineStyle', variant_lines{v});
    plot(ax(2), x, mean_hv(v,:), 'Color', variant_colors{v}, 'LineWidth', 2, 'LineStyle', variant_lines{v});
end

set(ax(1), 'YColor', [0.3, 0.3, 0.3], 'FontSize', 11);
set(ax(2), 'YColor', [0.3, 0.3, 0.3], 'FontSize', 11);
set(ax(1), 'XColor', [0.3, 0.3, 0.3], 'FontSize', 11);

ylabel(ax(1), 'Fitness', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(ax(2), 'Hypervolume', 'FontSize', 12, 'FontWeight', 'bold');
xlabel(ax(1), 'Iteration', 'FontSize', 12, 'FontWeight', 'bold');

title('Convergence Analysis: Proposed vs Ablation Variants', 'FontSize', 14, 'FontWeight', 'bold');

legend_handles = [];
legend_labels = [];
for v = 1:length(variants)
    legend_handles = [legend_handles, plot(nan, nan, 'Color', variant_colors{v}, 'LineWidth', 2, 'LineStyle', variant_lines{v})];
    legend_labels = [legend_labels, variant_labels{v}];
end
legend(legend_handles, legend_labels, 'Location', 'best', 'FontSize', 10);

grid(ax(1), 'on');
grid(ax(2), 'on');

break_point = 120;
annotation('arrow', [0.3, 0.45], [0.55, 0.7], 'Color', [1, 0, 0], 'LineWidth', 2);
annotation('textbox', [0.45, 0.72, 0.25, 0.08], 'String', 'Proposed escapes local optimum', ...
    'Color', [1, 0, 0], 'FontSize', 10, 'FontWeight', 'bold', ...
    'EdgeColor', [1, 0, 0], 'FaceAlpha', 0.1, 'BackgroundColor', [1, 1, 1]);

saveas(fig, fullfile(output_dir, 'convergence_hv.fig'));
saveas(fig, fullfile(output_dir, 'convergence_hv.png'));
fprintf('Convergence & HV plot saved to %s\n', output_dir);
close(fig);
end