function plot_boxplot(data_file, output_dir)
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

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4]);
set(gcf, 'Color', 'w');

variants = {'proposed', 'no_subpop', 'no_adaptive', 'no_levy', 'no_stop'};
variant_labels = {'Proposed', 'w/o Multi-Subpop', 'w/o Adaptive', 'w/o E-Levy', 'w/o MO-SmartStop'};

hv_data = cell(length(variants), 1);
for v = 1:length(variants)
    variant = variants{v};
    if isfield(results, variant)
        hv_data{v} = results.(variant).hv_values;
    else
        hv_data{v} = [];
    end
end

subplot(1, 2, 1);
boxplot(cell2mat(hv_data), 'Labels', variant_labels, 'Whisker', 1.5);
ylabel('Hypervolume', 'FontSize', 11, 'FontWeight', 'bold');
title('HV Distribution (30 Independent Runs)', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

h = findobj(gca, 'Tag', 'Box');
colors = lines(length(variants));
for i = 1:length(h)
    set(h(i), 'FaceColor', colors(i,:) * 0.7 + 0.3);
end

xtickangle(45);

subplot(1, 2, 2);
cov_data = cell(length(variants), 1);
for v = 1:length(variants)
    variant = variants{v};
    if isfield(results, variant)
        cov_data{v} = results.(variant).cov_high;
    else
        cov_data{v} = [];
    end
end

boxplot(cell2mat(cov_data), 'Labels', variant_labels, 'Whisker', 1.5);
ylabel('Coverage (%)', 'FontSize', 11, 'FontWeight', 'bold');
title('High-Priority Coverage Distribution', 'FontSize', 12, 'FontWeight', 'bold');
grid on;

h = findobj(gca, 'Tag', 'Box');
for i = 1:length(h)
    set(h(i), 'FaceColor', colors(i,:) * 0.7 + 0.3);
end

xtickangle(45);

saveas(fig, fullfile(output_dir, 'boxplot.fig'));
saveas(fig, fullfile(output_dir, 'boxplot.png'));
fprintf('Boxplot saved to %s\n', output_dir);
close(fig);
end