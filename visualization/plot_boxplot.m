function plot_boxplot(data_file, output_dir)
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

colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880
];

all_hv = [];
group_idx = [];

for i = 1:length(variants)
    v = variants{i};
    if isfield(results, v) && isfield(results.(v), 'hv_values')
        hv_data = results.(v).hv_values;
        all_hv = [all_hv; hv_data(:)];
        group_idx = [group_idx; i * ones(length(hv_data(:)), 1)];
    end
end

fig = figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

hold on;
for i = 1:length(variants)
    idx = (group_idx == i);
    if any(idx)
        b = boxchart(group_idx(idx), all_hv(idx));
        b.BoxFaceColor = colors(i, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end

xticks(1:5);
xticklabels(labels);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on;
ax = gca;
ax.GridLineStyle = '--';
ax.GridAlpha = 0.3;
box on;

saveas(fig, fullfile(output_dir, 'boxplot.fig'));
saveas(fig, fullfile(output_dir, 'boxplot.png'));
fprintf('Boxplot saved to %s\n', output_dir);
close(fig);
end