function table_str = plot_statistical_tests(data_file, output_dir)
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

variants = {'proposed', 'no_subpop', 'no_adaptive', 'no_levy', 'no_stop'};
variant_labels = {'Proposed', 'w/o Multi-Subpop', 'w/o Adaptive', 'w/o E-Levy', 'w/o MO-SmartStop'};

if ~isfield(results, 'proposed')
    fprintf('Error: proposed variant not found in data\n');
    table_str = '';
    return;
end

proposed_hv = results.proposed.hv_values;

n_variants = length(variants);
p_values = zeros(n_variants, 1);
h_stats = zeros(n_variants, 1);
mean_improvement = zeros(n_variants, 1);

for v = 1:n_variants
    variant = variants{v};
    if isfield(results, variant) && strcmp(variant, 'proposed')
        p_values(v) = 1.0;
        h_stats(v) = 0;
        mean_improvement(v) = 0;
    elseif isfield(results, variant)
        variant_hv = results.(variant).hv_values;
        [p_values(v), h_stats(v)] = ranksum(proposed_hv, variant_hv);
        mean_improvement(v) = (mean(proposed_hv) - mean(variant_hv)) / mean(variant_hv) * 100;
    else
        p_values(v) = NaN;
        h_stats(v) = NaN;
        mean_improvement(v) = NaN;
    end
end

fig = figure('Units', 'normalized', 'Position', [0.1, 0.1, 0.6, 0.4]);
set(gcf, 'Color', 'w');

ax = axes('Parent', fig);
axis off;

col_labels = {'Variant', 'Mean HV', 'Std HV', 'Improvement (%)', 'p-value', 'Significant?'};
table_data = cell(n_variants + 1, length(col_labels));
table_data{1, 1} = 'Proposed';
table_data{1, 2} = sprintf('%.4f', mean(proposed_hv));
table_data{1, 3} = sprintf('%.4f', std(proposed_hv));
table_data{1, 4} = '-';
table_data{1, 5} = '-';
table_data{1, 6} = 'N/A';

for v = 2:n_variants
    variant = variants{v};
    table_data{v, 1} = variant_labels{v};
    if isfield(results, variant)
        table_data{v, 2} = sprintf('%.4f', mean(results.(variant).hv_values));
        table_data{v, 3} = sprintf('%.4f', std(results.(variant).hv_values));
    else
        table_data{v, 2} = 'N/A';
        table_data{v, 3} = 'N/A';
    end
    table_data{v, 4} = sprintf('%.2f', mean_improvement(v));
    table_data{v, 5} = sprintf('%.4f', p_values(v));
    if h_stats(v) == 1
        table_data{v, 6} = 'Yes (*, p<0.05)';
    else
        table_data{v, 6} = 'No';
    end
end

column_format = {'left', 'center', 'center', 'center', 'center', 'left'};
column_width = {0.25, 0.12, 0.12, 0.15, 0.12, 0.2};

t = uitable('Parent', ax, ...
    'ColumnName', col_labels, ...
    'Data', table_data(2:end,:), ...
    'ColumnFormat', column_format, ...
    'ColumnWidth', column_width, ...
    'RowName', variant_labels(2:end), ...
    'FontSize', 10, ...
    'FontName', 'Times New Roman', ...
    'BackgroundColor', [1, 1, 1], ...
    'ForegroundColor', [0, 0, 0], ...
    'OuterPosition', [0.05, 0.1, 0.9, 0.85]);

title_str = 'Wilcoxon Signed-Rank Test: Proposed vs Ablation Variants (HV)';
title(title_str, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Times New Roman');

table_str = sprintf('Wilcoxon Test Results (Proposed HV: %.4f +/- %.4f)\n', mean(proposed_hv), std(proposed_hv));
table_str = [table_str, '\n'];
table_str = [table_str, sprintf('%-20s | %10s | %10s | %12s | %10s | %s\n', ...
    'Variant', 'Mean HV', 'Std HV', 'Improvement', 'p-value', 'Significant?')];
table_str = [table_str, repmat('-', 1, 85), '\n'];
for v = 2:n_variants
    table_str = [table_str, sprintf('%-20s | %10s | %10s | %11.2f%% | %10.4f | %s\n', ...
        variant_labels{v}, ...
        table_data{v, 2}, ...
        table_data{v, 3}, ...
        mean_improvement(v), ...
        p_values(v), ...
        table_data{v, 6})];
end

fprintf('\n%s\n', table_str);

saveas(fig, fullfile(output_dir, 'statistical_tests.fig'));
saveas(fig, fullfile(output_dir, 'statistical_tests.png'));
fprintf('Statistical tests saved to %s\n', output_dir);
close(fig);
end