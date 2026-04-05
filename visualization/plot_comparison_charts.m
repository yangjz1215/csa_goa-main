function plot_comparison_charts(data_file, output_dir)
if nargin < 1 || isempty(data_file)
    data_file = fullfile('..', 'experiments', 'comparison_results_para_Map1_Medium_*.mat');
    files = dir(data_file);
    if ~isempty(files)
        [~, idx] = sort([files.datenum], 'descend');
        data_file = fullfile(files(idx(1)).folder, files(idx(1)).name);
    else
        error('No comparison results file found');
    end
end
if nargin < 2 || isempty(output_dir)
    output_dir = fullfile('..', 'figures');
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

load(data_file);

fig1 = figure('Position', [100, 100, 900, 500]);
set(gcf, 'Color', 'w');

colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

algorithms = {'cSA_GOA', 'PSO', 'GA', 'GOA', 'cSA', 'GWO'};
labels = {'cSA-GOA (Proposed)', 'PSO', 'GA', 'GOA', 'cSA', 'GWO'};

subplot(1, 2, 1);
hold on;
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'convergence_curves')
        curves = results.(alg).convergence_curves;
        if iscell(curves) && ~isempty(curves)
            max_len = max(cellfun(@length, curves));
            avg_curve = zeros(1, max_len);
            n_runs = 0;
            for run = 1:min(length(curves), 30)
                if ~isempty(curves{run})
                    c = curves{run};
                    avg_curve(1:length(c)) = avg_curve(1:length(c)) + c;
                    n_runs = n_runs + 1;
                end
            end
            if n_runs > 0
                avg_curve = avg_curve / n_runs;
                if a == 1
                    lw = 2.5;
                else
                    lw = 1.8;
                end
                plot(1:length(avg_curve), avg_curve, '-', 'Color', colors(a,:), 'LineWidth', lw, 'DisplayName', labels{a});
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
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
box on;

subplot(1, 2, 2);
hold on;
hv_means = zeros(1, length(algorithms));
hv_stds = zeros(1, length(algorithms));
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'hv_values')
        vals = results.(alg).hv_values;
        hv_means(a) = mean(vals);
        hv_stds(a) = std(vals);
    end
end
errorbar(1:length(algorithms), hv_means, hv_stds, '-o', 'Color', [0.0000, 0.4470, 0.7410], ...
    'LineWidth', 2.0, 'MarkerSize', 8, 'MarkerFaceColor', [0.0000, 0.4470, 0.7410], 'Capsize', 6);
xlabel('Algorithm', 'FontWeight', 'bold');
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Comparison (Mean ± Std)', 'FontWeight', 'bold');
xticks(1:length(algorithms));
xticklabels(labels);
xtickangle(45);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 10, 'LineWidth', 1.2);
grid on;
ax = gca; ax.GridLineStyle = ':'; ax.GridAlpha = 0.5;
box on;

saveas(fig1, fullfile(output_dir, 'comparison_convergence.fig'));
saveas(fig1, fullfile(output_dir, 'comparison_convergence.png'));
close(fig1);

fig2 = figure('Position', [100, 100, 700, 550]);
set(gcf, 'Color', 'w');
hold on;

marker_styles = {'*', 'o', 's', '^', 'd', 'v'};
sizes = [120, 40, 40, 40, 40, 40];

for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'cov_high') && isfield(results.(alg), 'energies')
        cov_vals = results.(alg).cov_high;
        eng_vals = results.(alg).energies;
        if ~isempty(cov_vals) && ~isempty(eng_vals)
            if a == 1
                scatter(cov_vals, eng_vals/1000, sizes(a), marker_styles{a}, ...
                    'filled', 'MarkerFaceColor', colors(a,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
            else
                scatter(cov_vals, eng_vals/1000, sizes(a), marker_styles{a}, ...
                    'filled', 'MarkerFaceColor', colors(a,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);
            end
        end
    end
end

xlabel('Coverage (%)', 'FontWeight', 'bold');
ylabel('Energy Consumption (kJ)', 'FontWeight', 'bold');
legend(labels, 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 9);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on;
box on;

saveas(fig2, fullfile(output_dir, 'comparison_pareto.fig'));
saveas(fig2, fullfile(output_dir, 'comparison_pareto.png'));
close(fig2);

fig3 = figure('Position', [100, 100, 700, 500]);
set(gcf, 'Color', 'w');

hv_all = [];
group_idx = [];
for a = 1:length(algorithms)
    alg = algorithms{a};
    if isfield(results, alg) && isfield(results.(alg), 'hv_values')
        vals = results.(alg).hv_values;
        hv_all = [hv_all; vals(:)];
        group_idx = [group_idx; a * ones(length(vals(:)), 1)];
    end
end

box_colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880;
    0.6350, 0.0780, 0.1840
];

hold on;
for a = 1:length(algorithms)
    idx = (group_idx == a);
    if any(idx)
        b = boxchart(group_idx(idx), hv_all(idx));
        b.BoxFaceColor = box_colors(a, :);
        b.BoxFaceAlpha = 0.6;
        b.MarkerStyle = 'o';
        b.MarkerColor = [0.2, 0.2, 0.2];
        b.LineWidth = 1.5;
    end
end

xticks(1:length(algorithms));
xticklabels(labels);
ylabel('Hypervolume (HV)', 'FontWeight', 'bold');
title('HV Distribution (30 Independent Runs)', 'FontWeight', 'bold');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11, 'LineWidth', 1.2);
grid on;
ax = gca; ax.GridLineStyle = '--'; ax.GridAlpha = 0.3;
box on;

saveas(fig3, fullfile(output_dir, 'comparison_boxplot.fig'));
saveas(fig3, fullfile(output_dir, 'comparison_boxplot.png'));
close(fig3);

fprintf('Comparison charts saved to %s\n', output_dir);
fprintf('  - comparison_convergence.fig/png (Convergence + HV bar)\n');
fprintf('  - comparison_pareto.fig/png (Coverage vs Energy scatter)\n');
fprintf('  - comparison_boxplot.fig/png (HV robustness boxplot)\n');
end