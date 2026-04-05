function plot_pareto_comparison(data_file, output_dir)
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

fig = figure('Position', [100, 100, 700, 550]);
set(gcf, 'Color', 'w');
hold on;

variants_all = {'proposed', 'no_subpop', 'no_adaptive', 'no_levy', 'no_stop'};
labels_all = {'Proposed', 'w/o Subpop', 'w/o Adaptive', 'w/o E-Levy', 'w/o MO-Stop'};

all_colors = [
    0.8500, 0.3250, 0.0980;
    0.0000, 0.4470, 0.7410;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560;
    0.4660, 0.6740, 0.1880
];
all_markers = {'*', 'o', '^', 'd', 's'};
all_sizes = [120, 45, 45, 45, 45];

for i = 1:length(variants_all)
    variant = variants_all{i};
    if isfield(results, variant)
        r = results.(variant);
        if isfield(r, 'pareto_fronts') && ~isempty(r.pareto_fronts)
            pf_data = r.pareto_fronts;
            if iscell(pf_data) && length(pf_data) >= 1
                all_pf = [];
                for run = 1:min(length(pf_data), 10)
                    if ~isempty(pf_data{run})
                        all_pf = [all_pf; pf_data{run}];
                    end
                end
                if ~isempty(all_pf)
                    scatter(all_pf(:,1)*100, all_pf(:,2)/1000, all_sizes(i), all_markers{i}, ...
                        'filled', 'MarkerFaceColor', all_colors(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
                end
            end
        elseif isfield(r, 'cov_high') && isfield(r, 'energies')
            cov_vals = r.cov_high(:);
            eng_vals = r.energies(:);
            if ~isempty(cov_vals) && ~isempty(eng_vals)
                scatter(cov_vals*100, eng_vals/1000, all_sizes(i), all_markers{i}, ...
                    'filled', 'MarkerFaceColor', all_colors(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1.2);
            end
        end
    end
end

xlabel('Coverage (%)', 'FontWeight', 'bold');
ylabel('Energy Consumption (kJ)', 'FontWeight', 'bold');

legend(labels_all, 'Location', 'best', 'FontName', 'Times New Roman', 'FontSize', 10);

set(gca, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.2);
grid on;
box on;

saveas(fig, fullfile(output_dir, 'pareto_comparison.fig'));
saveas(fig, fullfile(output_dir, 'pareto_comparison.png'));
fprintf('Pareto comparison saved to %s\n', output_dir);
close(fig);
end