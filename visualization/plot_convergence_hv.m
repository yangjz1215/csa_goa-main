function plot_convergence_hv(data_file, output_dir)
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

fig = figure('Position', [100, 100, 800, 550]);
set(gcf, 'Color', 'w');

subplot(1, 2, 1);
hold on;
for v = 1:length(variants)
    variant = variants{v};
    if isfield(results, variant) && isfield(results.(variant), 'convergence_curves')
        curves = results.(variant).convergence_curves;
        if iscell(curves) && ~isempty(curves)
            max_len = max(cellfun(@length, curves));
            avg_curve = zeros(1, max_len);
            n_runs = 0;
            for run = 1:min(length(curves), 10)
                if ~isempty(curves{run})
                    c = curves{run};
                    avg_curve(1:length(c)) = avg_curve(1:length(c)) + c;
                    n_runs = n_runs + 1;
                end
            end
            if n_runs > 0
                avg_curve = avg_curve / n_runs;
                if v <= 5
                    lw = 2.0;
                else
                    lw = 2.5;
                end
                plot(1:length(avg_curve), avg_curve, '-', 'Color', colors{v, :}, 'LineWidth', lw, 'DisplayName', labels{v});
            end
        end
    end
end
xlabel('Generations', 'FontWeight', 'bold');
ylabel('Best Fitness', 'FontWeight', 'bold');
title('Convergence Curves', 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;
box on;

subplot(1, 2, 2);
hold on;
hv_means = zeros(1, length(variants));
hv_stds = zeros(1, length(variants));
for v = 1:length(variants)
    variant = variants{v};
    if isfield(results, variant) && isfield(results.(variant), 'hv_values')
        vals = results.(variant).hv_values;
        hv_means(v) = mean(vals);
        hv_stds(v) = std(vals);
    end
end
errorbar(1:5, hv_means, hv_stds, '-o', 'Color', [0.0000, 0.4470, 0.7410], 'LineWidth', 2.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.0000, 0.4470, 0.7410], 'Capsize', 6);
xlabel('Variant', 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Comparison (Mean ± Std)', 'FontWeight', 'bold');
xticks(1:5);
xticklabels(labels);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on;
ax = gca;
ax.GridLineStyle = ':';
ax.GridAlpha = 0.5;
box on;

saveas(fig, fullfile(output_dir, 'convergence_hv.fig'));
saveas(fig, fullfile(output_dir, 'convergence_hv.png'));
fprintf('Convergence & HV plot saved to %s\n', output_dir);
close(fig);
end